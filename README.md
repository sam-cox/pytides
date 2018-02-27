pytides
=======

## About

Pytides is small Python package for the analysis and prediction of tides. Pytides can be used to extrapolate the tidal behaviour at a given location from its previous behaviour. The method used is that of harmonic constituents, in particular as presented by P. Schureman in Special Publication 98. The fitting of amplitudes and phases is handled by Scipy's leastsq minimisation function. Pytides currently supports the constituents used by NOAA, with plans to add more constituent sets. It is therefore possible to use the amplitudes and phases published by NOAA directly, without the need to perform the analysis again (although there may be slight discrepancies for some constituents).

It is recommended that all interactions with pytides which require times to be specified are in the format of naive UTC datetime instances. In particular, note that pytides makes no adjustment for summertime or any other civil variations within timezones.

## Requirements

* Numpy
* Scipy

## Installation

```easy_install pytides```

or

```pip install pytides```

should do the trick.

Mainly for my own reference (sanity), to get pytides and its dependencies all working in a Debian (mint) virtualenv:
```
sudo apt-get install liblapack-dev libatlas-base-dev gfortran
export LAPACK=/usr/lib/liblapack.so
export ATLAS=/usr/lib/libatlas.so
export BLAS=/usr/lib/libblas.so
pip install numpy
pip install scipy
pip install pytides
```
and you'll probably want to
```
pip install matplotlib
```
although this won't install all the backends for matplotlib, which is a headache for another day ([this](http://www.stevenmaude.co.uk/2013/09/installing-matplotlib-in-virtualenv.html) looks promising).

## Usage

Pytides is in its infancy, and hasn't yet been fully documented. The best way to get started would be to read [this example](https://github.com/sam-cox/pytides/wiki/Example-Pytides-Usage).
After that, you might try [making your own tide table](https://github.com/sam-cox/pytides/wiki/How-to-make-your-own-Tide-Table-using-Python-and-Pytides), where you can also find a method for handling timezones.
You can find information about using NOAA published Harmonic Constituents directly [here](https://github.com/sam-cox/pytides/wiki/How-to-use-the-NOAA%27s-published-Harmonic-Constituents-in-Python-with-Pytides).

If you want to know *how* Pytides works, it would be best to read *P. Schureman, Special Publication 98*. Alternatively, there is [my attempt](https://github.com/sam-cox/pytides/wiki/Theory-of-the-Harmonic-Model-of-Tides) to explain it on the wiki (although it's a little mathematical and not yet complete).
It is certainly possible to use Pytides successfully without any knowledge of its methods.

## Contribution

I would welcome any help with Pytides. Particularly if you have knowledge of constituent data (including node factors) which other institutions/packages use. I can be reached at sam.cox@cantab.net or via github.

Please note however that Pytide's lack of attempts to infer constituents, or exclude similar constituents based on their synodic periods is an intentional design decision.
