=============================================================
Chromosome Location Annotator for HCX Hierarchies
=============================================================

This tool reads an HCX (hierarchical CX2) network, uses ``HCX::members`` and
``HCX::interactionNetworkUUID`` to pull genes from the referenced interaction
network, resolves chromosomes (human by default), tallies counts per
chromosome for every hierarchy node, adds those counts as node attributes, and
applies a pie-chart node style before returning an ``updateNetwork`` payload.

What it does
------------

* Reads HCX/CX2 input from a file or stdin
* Uses ``HCX::members`` + ``HCX::interactionNetworkUUID`` to fetch the referenced interaction network and pull genes from it
* Resolves genes to chromosomes using the bundled HGNC ``non_alt_loci_set.json`` map (or a custom map built via ``scripts/build_gene_chr_map.py``)
* Adds per-chromosome count attributes (``chr1_count`` … ``chr22_count``, ``chrX_count``, ``chrY_count``, ``chrM_count``)
* Adds summary attributes ``chromosomeCounts`` and ``chromosomeCountsJson``
* Applies a pie-chart node custom graphic driven by the chromosome count attributes
* Emits the modified CX2 wrapped in a single ``updateNetwork`` action

Dependencies
------------

* Python 3.11+
* ndex2

Install (for local testing)
---------------------------

.. code-block:: bash

   # from repository root
   python -m venv .venv
   source .venv/bin/activate
   pip install -r requirements.txt
   pip install -e .

Usage
-----

Read CX2 from a file and write updated network JSON to stdout:

.. code-block:: bash

   cytoscape-chromloc foo.cx2 --chrom-map path/to/non_alt_loci_set.json

Build a local gene→chromosome map (JSON/TSV) from HGNC:

.. code-block:: bash

   # download HGNC gene info from HGNC
   # grab json from total approved symbols, which includes chromosome info and is smaller than the full set of gene info
   # https://www.genenames.org/download/statistics-and-files/
   curl https://storage.googleapis.com/public-download-files/hgnc/json/json/non_alt_loci_set.json

Output format
-------------

The output JSON is a single action wrapping the modified network:

.. code-block:: json

   [{
     "action": "updateNetwork",
     "data": [ ... CX2 with chromosome attributes and pie style ... ]
   }]


Container build and usage
-------------------------

This app includes a simple Dockerfile under ``docker/`` so it can be run as a
CytoContainer-style service app.

Build image
~~~~~~~~~~~

From the repository root (same level as ``pyproject.toml``) run:

.. code-block:: bash

   cd docker

   # Build the image (same as the plain docker build command above)
   make build

Above will create `cytoscape-chromloc:latest` image. For additional configuration options
invoke `make` with no arguments

Basic usage
~~~~~~~~~~~

Show help for the service app inside the container:

.. code-block:: bash

   docker run --rm cytoscape-chromloc:latest

Run update on a local CX2 file (read-only bind mount assuming file is named foo.cx2):

.. code-block:: bash

   docker run --rm \
     -v "$(pwd)":/data \
     cytoscape-chromloc:latest \
     /data/hierarchy.cx2

The updated network JSON will be output to stdout


Integration with Cytocontainer
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Copy ``python_chromloc.json`` directory to the directory specified by 
``cytocontainer.algorithm.conf.dir`` in the configuration for your local installation of the
Cytoscape Container REST Service and restart the Cytoscape Container REST Service.
The app will be available at the endpoint ``/v1/algorithms/chromloc_annotator`` 

and can be invoked with a POST request containing a CX2 network in the body and the header ``Content-Type: application/json``.

The post should be JSON and look like this:

.. code-block:: json
 
  {
    "parameters": { 
                    "Updated By": "my tool"
     }, 
    "data": [{"CXVersion": "2.0", "hasFragments": false}, {"metaData": [{"elementCount": 1, "name": "attributeDeclarations"}, {"elementCount": 1, "name": "networkAttributes"}, {"elementCount": 1, "name": "nodes"}]}, {"attributeDeclarations": [{"networkAttributes": {"name": {"d": "string"}}, "nodes": {"name": {"d": "string"}, "represents": {"d": "string"}}}]}, {"networkAttributes": [{"name": "empty network"}]}, {"nodes": [{"id": 0, "x": 10, "y": 0, "v": {"name": "node 1", "represents": "representing node1"}}]}, {"edges": []}, {"status": [{"error": "", "success": true}]}]
  }

The response will be a JSON object containing id of the task

.. code-block:: json

  {
    "id": "some-uuid"
  }
