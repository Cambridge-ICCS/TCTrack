Machine Learning
================

The machine learning (ML) tracker detects tropical cyclones using a U-Net
segmentation model, and then links the detections into trajectories. The model
weights are downloaded from a HuggingFace Hub repository (or supplied locally)
rather than wrapping an external executable.

The model and its training pipeline are described in the reference `cyclone-track-ml
<https://github.com/robert-edwin-rouse/cyclone-track-ml>`_ repository.

For an overview of the functionalities, installation, and usage see the
:doc:`Machine Learning section of the documentation <../tracking-algorithms/machine_learning>`.

.. py:module:: tctrack.machine_learning

.. autoclass:: MLTracker
   :members:
   :undoc-members:
   :show-inheritance:
   :inherited-members:
   :exclude-members: model

.. autoclass:: MLParameters
   :members:
   :undoc-members:
   :show-inheritance:
   :inherited-members:

.. autoclass:: MLStitchParameters
   :members:
   :undoc-members:
   :show-inheritance:
   :inherited-members:

.. py:module:: tctrack.core.ml_tracker
   :no-index:

.. autoclass:: tctrack.core.ml_tracker.TCMLParameters
   :members:
   :undoc-members:
   :show-inheritance:

.. autoclass:: tctrack.core.ml_tracker.TCMLTracker
   :members:
   :undoc-members:
   :show-inheritance:
   :exclude-members: model
