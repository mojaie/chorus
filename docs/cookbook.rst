
Cookbook
====================================================================


Chemical structure drawing
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

SVG

::

    >>> from chorus.demo import MOL
    >>> from chorus import v2000reader as reader
    >>> from chorus.draw.svg import SVG
    >>> mol = reader.mol_from_text(MOL["demo"])
    >>> svg = SVG(mol)
    >>> svg.contents()
    '<svg xmlns="http://www.w3.org/2000/svg" xmlns:xlink="http://www.w3.org/1999/xlink" version="1.2" ...
    >>> svg.data_url_scheme()
    'data:image/svg+xml;base64,PHN2ZyB4bWxucz0iaHR0cDovL3d3dy53My5vcmcvMjAwMC9zdmciIHhtbG5zOnhsaW5rPS ...
    >>> svg.save("demo.svg")

.. image:: ../img/demo.svg
   :width: 300px


PNG (using matplotlib)

::

    >>> from chorus.demo import MOL
    >>> from chorus import v2000reader as reader
    >>> from chorus.draw.matplotlib import Matplotlib
    >>> mol = reader.mol_from_text(MOL["demo"])
    >>> mpl = Matplotlib(mol)
    >>> mpl.save("demo.png")

.. image:: ../img/demo.png
   :width: 300px



Chemical properties and descriptors
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

::

    >>> from chorus.demo import MOL
    >>> from chorus import molutil
    >>> from chorus import v2000reader as reader
    >>> mol = reader.mol_from_text(MOL["Paclitaxel"])
    >>> molutil.mw(mol)
    853.92
    >>> molutil.formula(mol)
    'C47H50NO14'
    >>> molutil.H_donor_count(mol)
    4
    >>> molutil.H_acceptor_count(mol)
    15
    >>> molutil.rotatable_count(mol)
    15
    >>> from chorus import wclogp
    >>> wclogp.wclogp(mol)
    3.53



* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
