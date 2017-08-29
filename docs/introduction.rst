
introduction
====================================================================

Chorus is a Python library for simple molecule modeling and cheminformatics analysis. Chorus now supports 2D molecular model based on graph theory and some descriptors which are used in the field of drug discovery.

Unlike RDKit, a front-runner project in this category, and other many integrated cheminformatics libraries, Chorus is designed to be simple and specialized for molecule itself hackability and portability.

Chorus molecule, atom and bond objects are pure Python objects which consist of dict-based object. It is easy to change the molecule component and properties for your experiment. This dict-based molecule model is easily converted to JSON like object which has high portablity to other platforms via HTTP access.

Chorus is not intended to be the alternative of other existing cheminformatics toolkits. Chorus highly depends on RDKit for its functionarity in cheminformatic analysis.



Features
----------

- Structure image export (PNG, SVG)
- Import from/export to .sdf, .mol
- Import from/export to RDKit molecule
- Molecular property calculation (MW, Chemical formula)
- Descriptors

  - H-bond donor/acceptor
  - Wildman-Crippen logP
  - Aromaticity

- Molecule graph topology (ring, scaffold, connectivity)
- Sub(super)structure search
- MCS with diameter restriction (MCS-DR) and graph-based local similarity (GLS)



Features (WIP)
-------------------------------

- Functional group descriptors
- Markush structue
- SMILES and 2D coordinate generation



Features (will never be implemented)
-------------------------------------

- Python 2 compatibility
- Fingerprint similarity
- And many of the features already available in RDKit



* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
