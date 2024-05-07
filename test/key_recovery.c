#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <time.h>
#include <math.h>

#include "../fips202.h"
#include "../sign.h"
#include "../packing.h"
#include "../polyvec.h"
#include "../polymatrix.h"
#include "../params.h"

#define MSGLEN ( 20 )

/*
 * The public key of EagleSign consists of (rho, A, E), where:
 * - rho: pseudorandom bitstring of SEEDBYTES bytes long,
 * - A: kxl matrix over R,
 * - E = (A + D) G^{-1} (mod q): kxl matrix over R.
 */
uint8_t public_key[CRYPTO_EAGLESIGN_PUBLICKEYBYTES], rho[SEEDBYTES];
polyvecl E[K], A[K];

/*
 * Given an unpacked signature (r, Z, W), return C = H(mu, r). The code is adapted from the
 * function `crypto_sign_verify()` in `sign.c` of the EagleSign3 reference implementation.
 */
void determine_c(polyvecl *C, uint8_t m[MSGLEN], uint8_t *r, polyvecl *Z, polyveck *W)
{
	uint8_t c[SEEDBYTES], mu[CRHBYTES];
	keccak_state state;
	uint16_t nonce_c = 0;

	/* Compute mu = CRH(H(pk), msg) */
	shake256(mu, SEEDBYTES, public_key, CRYPTO_EAGLESIGN_PUBLICKEYBYTES);
	shake256_init(&state);
	shake256_absorb(&state, mu, SEEDBYTES);
	shake256_absorb(&state, m, MSGLEN);
	shake256_finalize(&state);
	shake256_squeeze(mu, CRHBYTES, &state);

	/* Call the random oracle and Compute C = H(mu,r) dans B_tau^l*/
	shake256_init(&state);
	shake256_absorb(&state, r, SEEDBYTES);
	shake256_absorb(&state, mu, CRHBYTES);
	shake256_finalize(&state);
	shake256_squeeze(c, SEEDBYTES, &state);
	polyvecl_challenge_y1_c(C, c, &nonce_c, 1);

	/* Put C in coefficient form. */
	polyvecl_invntt_tomont(C);
}

/*
 * Given a function that can generate random (message, signature)-pairs,
 * recover the matrices G, and D, which amounts to secret key recovery.
 */
void private_key_recovery(
	polyvecl pred_G[L], polyvecl pred_D[K], int num_sigs,
	void (*gensig)(uint8_t[MSGLEN], uint8_t[CRYPTO_EAGLESIGN_BYTES]))
{
	// message (used to call `gensig`), and raw signature
	uint8_t msg[MSGLEN], sig[CRYPTO_EAGLESIGN_BYTES];

	// Unpacked signature
	uint8_t sig_r[SEEDBYTES];
	polyvecl sig_Z, sig_C;
	polyveck sig_W;

	// Statistics on the signatures
	long long num_traces[L] = {}, traceG[L][L][N] = {}, traceD[K][L][N] = {};
	double recG[N], recD[N];

	// Expand the public key
	unpack_pk(rho, E, public_key);
	polyvec_matrix_expand(A, rho);

	for (int i = 0; i < num_sigs; i++) {
		gensig(msg, sig);
		assert(unpack_sig(sig_r, &sig_Z, &sig_W, sig) == 0);
		assert(!polyvec_chknorms(&sig_Z, &sig_W));
		determine_c(&sig_C, msg, sig_r, &sig_Z, &sig_W);

		for (int lu = 0; lu < L; lu++) {

			for (int idx_c = 0; idx_c < N; idx_c++) {
				if (sig_C.vec[lu].coeffs[idx_c] == 0) {
					continue;
				}

				long long c0 = sig_C.vec[lu].coeffs[idx_c];
				num_traces[lu]++;

				// Gather statistics on G
				for (int lz = 0; lz < L; lz++) {
					for (int j = 0; j < N; j++) {
						int idx_z = idx_c + j, sign = 1;
						if (idx_z >= N) idx_z -= N, sign = -sign;
						long long zj = sign * (long long)sig_Z.vec[lz].coeffs[idx_z];
						traceG[lz][lu][j] += c0 * zj;
					}
				}

				// Gather statistics on D
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

	// Recover the secret key
	for (int lu = 0; lu < L; lu++) {
		// Recover G
		for (int lz = 0; lz < L; lz++) {
			double sq_norm = 0;
			for (int i = 0; i < N; i++) {
				recG[i] = ((double)traceG[lz][lu][i]) / num_traces[lu];
				sq_norm += recG[i] * recG[i];
			}
			double renorm = sqrt(2 * N / (3 * sq_norm));
			for (int i = 0; i < N; i++) {
				pred_G[lz].vec[lu].coeffs[i] = lround(recG[i] * renorm);
			}
		}

		// Recover D
		for (int lw = 0; lw < K; lw++) {
			double sq_norm = 0;
			for (int i = 0; i < N; i++) {
				recD[i] = ((double)traceD[lw][lu][i]) / num_traces[lu];
				sq_norm += recD[i] * recD[i];
			}
			double renorm = sqrt(2 * N / (3 * sq_norm));
			for (int i = 0; i < N; i++) {
				pred_D[lw].vec[lu].coeffs[i] = lround(recD[i] * renorm);
			}
		}
	}
}

/*
 * The private key of EagleSign consists of (rho, A, G, D), where:
 * - rho, A: similar to the public key,
 * - G: lxl matrix over R,
 * - D: kxl matrix over R.
 */
uint8_t private_key[CRYPTO_EAGLESIGN_SECRETKEYBYTES];

void generate_msg_sig(uint8_t msg[MSGLEN], uint8_t sig[CRYPTO_EAGLESIGN_BYTES]) {
	size_t siglen;

	for (int i = 0; i < MSGLEN; i++) {
		msg[i] = rand() & 0xFF;
	}

	assert(crypto_sign_signature(sig, &siglen, msg, MSGLEN, private_key) == 0);
	assert(siglen == CRYPTO_EAGLESIGN_BYTES);
}

/*
 * Check whether the G and D matrices correctly recovered.
 */
void check_key_recovery(const polyvecl pred_G[L], const polyvecl pred_D[K],
	int *num_wrong_G, int *num_wrong_D)
{
	uint8_t rho[SEEDBYTES], tr[SEEDBYTES];
	polyvecl G[L], D[K];
	unpack_sk(rho, tr, G, D, private_key);

	*num_wrong_G = *num_wrong_D = 0;
	for (int lu = 0; lu < L; lu++) {
		for (int i = 0; i < N; i++) {
			for (int lz = 0; lz < L; lz++) {
				*num_wrong_G += G[lz].vec[lu].coeffs[i] != pred_G[lz].vec[lu].coeffs[i];
			}
			for (int lw = 0; lw < K; lw++) {
				*num_wrong_D += D[lw].vec[lu].coeffs[i] != pred_D[lw].vec[lu].coeffs[i];
			}
		}
	}
}

int main(int argc, char **argv) {
	polyvecl pred_G[L] = {}, pred_D[K] = {};
	int num_sigs = argc > 1 ? atoi(argv[1]) : 1000, num_wrong_G, num_wrong_D;
	srand(time(NULL));

	if (argc > 2) {
		int num_runs = atoi(argv[2]), success_G = 0, success_D = 0;
		for (int run = 0; run < num_runs; run++) {
			crypto_sign_keypair(public_key, private_key);
			private_key_recovery(pred_G, pred_D, num_sigs, generate_msg_sig);

			// Check whether recovery was successful.
			check_key_recovery(pred_G, pred_D, &num_wrong_G, &num_wrong_D);
			success_G += num_wrong_G == 0;
			success_D += num_wrong_D == 0;
		}

		printf("Success rate recovering G: %.3f%%\n", 100.0 * success_G / num_runs);
		printf("Success rate recovering D: %.3f%%\n", 100.0 * success_D / num_runs);
	} else {
		// Key Generation
		crypto_sign_keypair(public_key, private_key);

		// Recover G
		private_key_recovery(pred_G, pred_D, num_sigs, generate_msg_sig);

		// Check whether recovery was successful.
		check_key_recovery(pred_G, pred_D, &num_wrong_G, &num_wrong_D);
		if (num_wrong_G == 0 && num_wrong_D == 0) {
			printf("Key recovery attack was successful.\n");
		} else {
			printf("%4d/%4d coefficients of G are incorrectly recovered.\n", num_wrong_G, N * L * L);
			printf("%4d/%4d coefficients of D are incorrectly recovered.\n", num_wrong_D, N * L * K);
		}
	}
	return 0;
}
