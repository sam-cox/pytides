pytides
=======
####About
Pytides is small Python package for the analysis and prediction of tides. Pytides can be used to extrapolate the tidal behaviour at a given location from its previous behaviour. The method used is that of harmonic constituents, in particular as presented by P. Schureman in Special Publication 98. The fitting of amplitudes and phases is handled by Scipy's leastsq minimisation function. Pytides currently supports the constituents used by NOAA, with plans to add more constituent sets. It is therefore possible to use the amplitudes and phases published by NOAA directly, without the need to perform the analysis again (although there may be slight discrepancies for some constituents).

It is recommended that all interactions with pytides which require times to be specified are in the format of naive UTC datetime instances. In particular, note that pytides makes no adjustment for summertime or any other civil variations within timezones.



####Installation

```easy_install pytides```

or

```pip install pytides```

should do the trick.

####Usage
Pytides is in its infancy, and hasn't yet been properly documented. The best way to get started would be to read [this example](https://github.com/sam-cox/pytides/wiki/Example-Pytides-Usage).

####Contribution
I would welcome any help with Pytides. Particularly if you have knowledge of constituent data (including node factors) which other institutions/packages use. I can be reached at sam.cox@cantab.net or via github.

Please note however that Pytide's lack of attempts to infer constituents, or exclude similar constituents based on their synodic periods is an intentional design decision.
