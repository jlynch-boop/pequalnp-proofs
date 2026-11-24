# Partial results for L-gct-210

L-gct-210 asks whether a complex linear space of $5\times5$ matrices with $\operatorname{perm}(J+T)=120$ for every direction $T$ has dimension at most eight. Here $J$ is the all-ones matrix. The unrestricted conjecture remains open. The bounds below are mathematical proofs for the stated families, with exact arithmetic checks of selected identities. [Project desk](https://data.pequalnp.com/data/desk/gaps.json)

## Dimension bounds

- **Common zero block.** The dimension is at most eight when every direction vanishes on the same $2\times3$ block. The proof reduces the remaining case to five-dimensional spaces of $3\times3$ matrices of permanent zero and uses the cofactors of the last two columns to finish. [Proof](zero-block-bound.txt) · [Exact identity checks](check_zero_block_identities.py) · [Chan and Ilten, Section 7B](https://msp.org/ant/2015/9-3/ant-v9-n3-p05-s.pdf)
- **Full variation in one column or row.** The dimension is at most eight when the space contains all sum-zero directions supported in one specified column or row. In the column case, varying that column makes all five cofactors equal to 24. The projected space then has dimension at most four, by the quadratic identities and a cubic and quartic calculation on the kernel of its row-sum map. [Proof](full-column-variance-bound.txt)
- **Standard example and centered core.** The standard seven-dimensional example cannot be enlarged by any linear direction while retaining the constant permanent. The bound of eight also holds in the 15-dimensional family with a doubly centered $4\times4$ core, sum-zero borders, and zero corner. [Standard example](standard-example.txt) · [Centered core and general constraints](arbitrary-spaces.txt)

## Fourier search

We checked coordinate spans of Fourier modes from every pair of row and column cycle types, up to transposition. Across 28 pairs, the largest passing span modulo 61 had dimension six. Mapping the cyclotomic coefficients into this field preserves any exact identity, so each failed modular test rules out the corresponding complex span. [Search code](search_fourier_spans.py) · [Recorded output](fourier-scan.txt)

Run the checks from this directory with Python 3. The Fourier search also uses NumPy.

```sh
python3 verify_star.py
python3 check_zero_block_identities.py
python3 search_fourier_spans.py
```

## Remaining column case

Three independent directions supported in one column constrain the cofactor vector to a two-dimensional space. When the five pencil parameters are distinct and the row-sum map is surjective, the remaining case reduces to a rank-three symmetric pencil on a four-dimensional space. The cubic equations exclude one explicit realization. The other realizations and repeated-parameter cases remain open. [Reduction and example](three-dimensional-column-pencil.txt)
