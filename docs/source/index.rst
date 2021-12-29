.. GoValve-Root documentation master file, created by
   sphinx-quickstart on Thu Nov  4 16:08:12 2021.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to GoValve-Root's documentation!
========================================

**GoValve-Root** is an ITK-SNAP distributed segmentation service (DSS) that computes aortic root
strain over the cardiac cycle from a 4D image series.

.. toctree::
   :maxdepth: 2
   :caption: Contents:

Prerequisites
=============
<<<<<<< HEAD
* `ITK-SNAP <itksnap.org>` version 4.0 or later
Note that earlier versions of ITK-SNAP do not handle 4D segmentations.
=======
*`ITK-SNAP <itksnap.org>` 4.0 or later

>>>>>>> ccd2497a1813488e87bda41ae439571f3f370ce4

DSS Pipeline Overview
=====================
As with any DSS implemented in ITK-SNAP, the GoValve-Root service involves three layers of communication:

**client**
  The user supplies the input to the algorithm by loading a 4D image of the aortic root in the ITK-SNAP GUI (version 4.0 or later), creating a 3D segmentation of the aortic root in a "reference frame" (a single 3D image volume in the 4D series), and tagging relevant time points in the cardiac cycle.

**middleware**
  The ITK-SNAP workspace created by the client is submitted to the ITK-SNAP DSS middleware layer, which orchestrates communication between the client and the service provider that carries out the image analysis algorithm. The main production DSS middleware runs at https://dss.itksnap.org. The user must create an account and sign in to dss.itksnap.org before submitting an image for segmentation and strain analysis.

**service**
  The server that runs the strain analysis algorithm receives the image from the DSS middleware layer and returns the output 4D aortic root segmentation, mesh series, and strain information.  

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
