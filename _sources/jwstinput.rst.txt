Possible Instrument Input Params
================================

Instrument dictionaries use Pandeia's lowercase configuration values. Load a
template with ``jdi.load_mode_dict('NIRSpec G395H')``, then edit its
``configuration`` keys.
The values below are the combinations exposed by PandExo's web interface;
confirm final timing and detector settings in the current APT.

NIRSpec 
-------

Disperser/filter combinations
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

- ``g140m`` or ``g140h`` with ``f070lp`` or ``f100lp``
- ``g235m`` or ``g235h`` with ``f170lp``
- ``g395m`` or ``g395h`` with ``f290lp``
- ``prism`` with ``clear``

All grating modes support ``sub2048``, ``sub1024a``, ``sub1024b``, and
``sub512``. PRISM/CLEAR does not support ``sub1024a``; it additionally supports
the multistripe subarrays ``s256m2_prm``, ``s128m4_prm``, ``s64m8_prm``, and
``s32m16_prm``. Multistripe subarrays are valid only with PRISM/CLEAR.

NIRISS
------

NIRISS SOSS supports ``substrip256``, ``substrip96``,
``sub17stripe_soss``, ``sub60stripe_soss``, ``sub204stripe_soss``, and
``sub680stripe_soss``. Order 1 is the default; order 2 extraction is available
with ``substrip256`` by setting ``inst_dict['strategy']['order'] = 2``.

NIRCam
------

Long-wave grism time series use ``grismr`` with ``f322w2`` or ``f444w`` and
``subgrism64``, ``subgrism128``, or ``subgrism256``. One- and four-output
variants are available.

The ``NIRCam DHS`` template pairs one short-wave filter (``f070w``, ``f090w``,
``f115w``, ``f150w2``, or ``f200w``) with ``f322w2`` or ``f444w``. Supported
DHS subarrays are ``sub41s1_2-spectra``, ``sub82s2_4-spectra``,
``sub164s4_8-spectra``, and ``sub260s4_8-spectra``. Short- and long-wave
channels are observed simultaneously but returned one channel per calculation.

MIRI
----

MIRI slitless LRS uses ``p750l`` and ``slitlessprism``,
``slitlessprism_ip``, or ``slitlessprism_ips``. MIRI LRS slit mode can be
simulated with ``subslit`` or ``full``, but it is not currently offered for
JWST time-series observations and PandExo reports a warning.

Groups and readouts
-------------------

Set ``inst_dict['configuration']['detector']['ngroup'] = 'optimize'`` to let
PandExo select a group count, or enter an integer valid for the selected
Pandeia detector configuration. NIRCam and DHS templates also accept
``readout_pattern = 'optimize'`` to account for saturation and estimated data
excess when choosing a readout.

Use ``jdi.subarrays()``, ``jdi.dispersers()``, and ``jdi.filters()`` to inspect
the values known to the installed PandExo release.
