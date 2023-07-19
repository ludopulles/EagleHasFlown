# EagleHasFlown
A key recovery attack on EagleSign.

## Basic explanation

In EagleSign, signatures contain `z = Gu (mod q)`, with `u = y1 + c`, where `G, y1` and `c` are ternary.
Both `z` and `c` can be obtained from a (message, signature) pair, and by [Mehdi Tibouchi](https://groups.google.com/a/list.nist.gov/g/pqc-forum/c/zas5PLiBe6A/m/ik5IpWd7BQAJ)'s observation, this leaks information about the matrix G.

For the simple case, suppose `L = K = 1` so we do not consider matrices yet, but note that we still work over the number ring `Z[X]/(X^n+1)`.
Consider many signatures where `c` has first coefficient equal to `1`.
In that case, `z = G*y1 + G*c` will have an average value of `G * 1`, since `y1` is zero on average, and similarly all the other coefficients of `c` are also zero.
Given many signatures with a `c` that has a `1` in its first coefficient, we can estimate the average value of `z`, thus revealing the value of `G`.
However, we can improve the situation by also considering `c` with first coefficient equal to `-1`, as in that case `z` is on average `-G`, or equivalently, `-z` is on average `G`.

Further improvement can be done by considering other coefficients of `c`.
Namely, write `c = c_0 + c_1 X + ... + c_{n-1} X^{n-1}`.
By selecting all signatures with `c_k` nonzero, note that `G*(y1 + c)` is on average `G * (c_k X^k)`, or equivalently, `z c_k X^{-k}` is on average `G c_k^2 = G` (using c_k = -1,1).

Now given enough signatures under arbitrary messages, we succesfully recover `G`.
In particular, already 300 signatures seems to be enough to recover G when using EagleSign3, while 500 signatures are enough to confidently recover G when using EagleSign5.

## Matrix case

This attack is easily extended to the matrix situation (`L > 1` or `K > 1`), because the expectation values of `y1` and `c` are still zero vectors.
Now to recover the `(i,j)`th entry `G_{i,j}` from the LxL matrix `G`, observe that `z_i = \sum_{j=1}^L G_{i,j} (y1_j + c_j)` so we can use the above with now `z_i` and `c_j` in the role of `z` and `c`.

## Recovering D

The secret KxL matrix D can recovered with the equation `w = y2 - D u (mod q)`, where y2 is a vector with small (abs. value <= 64) coefficients.
Again, `w` and `c = u - y1` can be obtained from a (message, signature) pair.
Taking the minus sign here into account, we can still basically use the above method to recover D, because the y2 vector has coefficients much smaller than q, so there is still no wrapping around.
However, now roughly ~2000 signatures is enough to confidently recover D.

## Usage

To run the attack:

- `make`
- Execute `test/key_recovery3 [num_signatures]` for key recovery on EagleSign3
- Execute `test/key_recovery5 [num_signatures]` for key recovery on EagleSign5

