---
jupytext:
  formats: ipynb,md:myst
  text_representation:
    extension: .md
    format_name: myst
    format_version: 0.13
    jupytext_version: 1.13.8
kernelspec:
  display_name: Python 3 (ipykernel)
  language: python
  name: python3
---

Run the tests in the root directory:

    fpm test

This will generate the `data_*.txt` files that we can then plot below.

```{code-cell} ipython3
%pylab inline
```

```{code-cell} ipython3
D = loadtxt("../data_harmonic_schroed.txt")
x = D[0,:]
#n = size(D,0)
n = 8
figure(figsize=(20,12))
for i in range(1,n):
    plot(x, D[i,:], "-", label=f"{i}")
xlim([0,5])
xlabel("r [a.u.]")
ylabel("wavefunction [a.u.]")
title("Schrödinger Harmonic Oscillator Wavefunctions")
legend()
savefig("harmonic_schroed.pdf")
show()
```

```{code-cell} ipython3
D = loadtxt("../data_harmonic_dirac.txt")
x = D[0,:]
#n = size(D,0)
n = 8
figure(figsize=(20,12))
for i in range(1,n):
    plot(x, D[i,:], "-", label=f"{i}")
xlim([0,5])
xlabel("r [a.u.]")
ylabel("wavefunction [a.u.]")
title("Dirac Harmonic Oscillator Wavefunctions")
legend()
savefig("harmonic_dirac.pdf")
show()
```

```{code-cell} ipython3
D = loadtxt("../data_coulomb_schroed.txt")
x = D[0,:]
#n = size(D,0)
n = 8
figure(figsize=(20,12))
for i in range(1,n):
    semilogx(x, D[i,:]*x, "-", label=f"{i}")
xlim([None,1.5])
xlabel("r [a.u.]")
ylabel("wavefunction [a.u.]")
title("Schrödinger Coulomb Wavefunctions")
legend()
savefig("coulomb_schroed.pdf")
show()
```

```{code-cell} ipython3
D = loadtxt("../data_coulomb_dirac.txt")
x = D[0,:]
#n = size(D,0)
n = 8
figure(figsize=(20,12))
for i in range(1,n):
    semilogx(x, D[i,:]*x, "-", label=f"{i}")
xlim([None,1.5])
xlabel("r [a.u.]")
ylabel("wavefunction [a.u.]")
title("Dirac Coulomb Wavefunctions")
legend()
savefig("coulomb_dirac.pdf")
show()
```
