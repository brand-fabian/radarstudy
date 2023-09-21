.. pybwa documentation master file, created by
   sphinx-quickstart on Wed May 26 16:05:59 2021.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to cybwa's documentation!
=================================

.. toctree::
   :maxdepth: 2
   :caption: Contents:
   
.. automodule:: bwa
   :members:
   :undoc-members:
   :show-inheritance:
   :exclude-members: SequenceListForwardIterator, SequenceListView, SequenceList
   
   .. autoclass:: SequenceList
     :members:
     :undoc-members:
     :special-members: __iter__, __getitem__
     
   .. autoclass:: Bwa
     :members:
     :undoc-members:
   
   .. autoclass:: MEM_F
     :members:
     :undoc-members:
   
     Bwamem Option flags.

     Use these option flags to signal specific conditions to the bwa mem
     algorithm. Can be set using :func:`~bwa.BwamemOptions.set_flag`.

   .. autoclass:: BWTALGO
      :members:
      :undoc-members:
      
      Index generation algorithms.
      
   .. autoclass:: BWA_IDX
      :members:
      :undoc-members:
      
      BWA index load algorithm/type.
      
   .. autoclass:: SequenceListForwardIterator
      :members:
      :undoc-members:
      :special-members: __next__
      
   .. autoclass:: SequenceListView
      :members:
      :undoc-members:
      :special-members: __iter__, __getitem__


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
