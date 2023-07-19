#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <time.h>
#include <math.h>

#include "../fips202.h"
#include "../sign.h"
#include "../packing.h"
#include "../poly.h"
#include "../polyvec.h"
#include "../polymatrix.h"
#include "../params.h"
#include "../api.h"

#define MSGLEN ( 20 )

void print_G_and_D(const uint8_t *private_key)
{
	uint8_t rho[SEEDBYTES], tr[SEEDBYTES];
	polyvecl G[L], D[K];
	unpack_sk(rho, tr, G, D, private_key);
 
	// EagleSign3: L = K = 1

	// ================================================================================
	// Print G
	// ================================================================================
	for (int i = 0; i < L; i++) {
		for (int j = 0; j < L; j++) {
			printf("G[%d][%d] = (", i, j);
			for (int k = 0; k < N; k++) {
				if (k % 64 == 0) printf("\n");
				printf("%2d", G[i].vec[j].coeffs[k]);
				if (k != N-1) printf(",");
			}
			printf(")\n");
		}
	}

	return;

	// ================================================================================
	// Print D
	// ================================================================================
	for (int i = 0; i < K; i++) {
		for (int j = 0; j < L; j++) {
			printf("D[%d][%d] = (", i, j);
			for (int k = 0; k < N; k++) {
				if (k % 64 == 0) printf("\n");
				printf("%2d", D[i].vec[j].coeffs[k]);
				if (k != N-1) printf(",");
			}
			printf(")\n");
		}
	}
}

// Public key:
uint8_t public_key[CRYPTO_EAGLESIGN_PUBLICKEYBYTES], rho[SEEDBYTES];
polyvecl E[K];
polyvecl A[K]; // matrix A is even more 'unpacked'

void determine_c(polyvecl *C, uint8_t *m, size_t mlen, uint8_t *r, polyvecl *Z, polyveck *W)
{
	uint8_t c[SEEDBYTES], mu[CRHBYTES];
	keccak_state state;
	uint16_t nonce_c = 0;

/*
	unpack_pk(rho, E, pk);

	unpack_sig(r, &Z, &W, sig);

	if (polyvec_chknorms(Z, W))
		return -1;
*/

	/* Applying NTT Transformation*/
	polymatrix_ntt_k_l(E);
	polyvecl_ntt(Z);
	polyveck_ntt(W);

	/* Compute mu = CRH(H(pk), msg) */
	shake256(mu, SEEDBYTES, public_key, CRYPTO_EAGLESIGN_PUBLICKEYBYTES);
	shake256_init(&state);
	shake256_absorb(&state, mu, SEEDBYTES);
	shake256_absorb(&state, m, mlen);
	shake256_finalize(&state);
	shake256_squeeze(mu, CRHBYTES, &state);

	/* Call the random oracle and Compute C = H(mu,r) dans B_tau^l*/
	shake256_init(&state);
	shake256_absorb(&state, r, SEEDBYTES);
	shake256_absorb(&state, mu, CRHBYTES);
	shake256_finalize(&state);
	shake256_squeeze(c, SEEDBYTES, &state);
	polyvecl_challenge_y1_c(C, c, &nonce_c, 1);

	polyvecl_invntt_tomont(C);
}

void private_key_recovery(polyvecl *predicted_G, polyvecl *predicted_D, int num_sigs, void (*gensig)(uint8_t*, size_t*, uint8_t*))
{
	// ================================================================================
	uint8_t msg[MSGLEN], sig[CRYPTO_EAGLESIGN_BYTES];
	size_t msglen;

	// Unpacked signature:
	uint8_t sig_r[SEEDBYTES];
	polyvecl sig_Z, sig_C;
	polyveck sig_W;

	// Statistics on the signatures:
	long long num_traces[L] = {}, traceG[L][L][N] = {}, traceD[K][L][N] = {};
	double recG[N], recD[N];


	// ================================================================================
	// It's quite likely that you can extend the attack to L = 2.
	// For now, the goal is only to illustrate the weakness for L = 1.
	unpack_pk(rho, E, public_key);
	/* Expand matrix A in NTT form*/
	polyvec_matrix_expand(A, rho);

	for (int i = 0; i < num_sigs; i++) {
		gensig(msg, &msglen, sig);

		assert(unpack_sig(sig_r, &sig_Z, &sig_W, sig) == 0);
		assert(!polyvec_chknorms(&sig_Z, &sig_W));

		determine_c(&sig_C, msg, msglen, sig_r, &sig_Z, &sig_W);
		polyvecl_invntt_tomont(&sig_Z);
		polyveck_invntt_tomont(&sig_W);

		for (int lu = 0; lu < L; lu++) {

			for (int idx_c = 0; idx_c < N; idx_c++) {
				if (sig_C.vec[lu].coeffs[idx_c] != 0) {
					long long c0 = sig_C.vec[lu].coeffs[idx_c];
					num_traces[lu]++;

					// Gather statistics on G:
					for (int lz = 0; lz < L; lz++) {
						for (int j = 0; j < N; j++) {
							int idx_z = idx_c + j, sign = 1;
							if (idx_z >= N) idx_z -= N, sign = -sign;
							long long zj = sign * (long long)sig_Z.vec[lz].coeffs[idx_z];
							traceG[lz][lu][j] += c0 * zj;
						}
					}

					// Gather statistics on D:
					for (int lw = 0; lw < K; lw++) {
						for (int j = 0; j < N; j++) {
							int idx_w = idx_c + j, sign = 1;
							if (idx_w >= N) idx_w -= N, sign = -sign;
							long long wj = -sign * (long long)sig_W.vec[lw].coeffs[idx_w];
							traceD[lw][lu][j] += c0 * wj;
						}
					}

				}
			}
		}
	}

	// Recover the secret key.
	for (int lu = 0; lu < L; lu++) assert(num_traces[lu] == TAU * num_sigs);
	// printf("Took %lld traces\n", num_traces[0]);

	for (int lu = 0; lu < L; lu++) {

		// Recover G:
		for (int lz = 0; lz < L; lz++) {
			double norm_recG = 0;
			for (int i = 0; i < N; i++) {
				recG[i] = ((double)traceG[lz][lu][i]) / num_traces[lu];
				norm_recG += recG[i] * recG[i];
			}
			for (int i = 0; i < N; i++) {
				// Formula used from Mehdi's code:
				recG[i] *= sqrt(2 * N / (3 * norm_recG));
				predicted_G[lz].vec[lu].coeffs[i] = lround(recG[i]);
			}
		}

		// Recover D:
		for (int lw = 0; lw < K; lw++) {
			double norm_recD = 0;
			for (int i = 0; i < N; i++) {
				recD[i] = ((double)traceD[lw][lu][i]) / num_traces[lu];
				norm_recD += recD[i] * recD[i];
			}
			for (int i = 0; i < N; i++) {
				// Formula used from Mehdi's code:
				recD[i] *= sqrt(2 * N / (3 * norm_recD));
				predicted_D[lw].vec[lu].coeffs[i] = lround(recD[i]);
			}
		}

	}
}

void compare_keys(const uint8_t *private_key, const polyvecl *predicted_G, const polyvecl *predicted_D)
{
	uint8_t rho[SEEDBYTES], tr[SEEDBYTES];
	polyvecl G[L], D[K];
	unpack_sk(rho, tr, G, D, private_key);
 
	for (int lu = 0; lu < L; lu++) {

		// Show results for G:
		for (int lz = 0; lz < L; lz++) {
			int n_correct = 0;
			for (int i = 0; i < N; i++) {
				int aGi = G[lz].vec[lu].coeffs[i];
				int pGi = predicted_G[lz].vec[lu].coeffs[i];
				n_correct += (aGi == pGi);
			}
			printf("%d out of %d coefficients of G[%d,%d] are recovered succesfully.\n", n_correct, N, lz, lu);
		}

		// Show results for D:
		for (int lw = 0; lw < L; lw++) {
			int n_correct = 0;
			for (int i = 0; i < N; i++) {
				int aDi = D[lw].vec[lu].coeffs[i];
				int pDi = predicted_D[lw].vec[lu].coeffs[i];
				n_correct += (aDi == pDi);
			}
			printf("%d out of %d coefficients of D[%d,%d] are recovered succesfully.\n", n_correct, N, lw, lu);
		}
	}

}



uint8_t private_key[CRYPTO_EAGLESIGN_SECRETKEYBYTES];

void generate_msg_sig(uint8_t *msg, size_t *msglen, uint8_t sig[CRYPTO_EAGLESIGN_BYTES]) {
	size_t siglen;

	*msglen = MSGLEN;
	for (int i = 0; i < MSGLEN; i++) {
		msg[i] = rand() & 0xFF;
	}

	assert(crypto_sign_signature(sig, &siglen, msg, MSGLEN, private_key) == 0);
	assert(siglen == CRYPTO_EAGLESIGN_BYTES);
}

signed main(int argc, char **argv) {
	int num_sigs = argc > 1 ? atoi(argv[1]) : 1000;
	polyvecl predicted_G[L], predicted_D[K];

	srand(time(NULL));

	// Key Generation:
	crypto_sign_keypair(public_key, private_key);

	// Recover G:
	private_key_recovery(predicted_G, predicted_D, num_sigs, generate_msg_sig);

	// Print the actual private key.
	// print_G_and_D(private_key);

	// Compare G and predicted_G:
	compare_keys(private_key, predicted_G, predicted_D);

	return 0;
}
