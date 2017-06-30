Output Data
===========

The output data is create in a strictly structured manor.
All paths are given relative to the base_dir parameter, the next level gives the major.minor version of the code that created the data and the next level the pipeline name.

In luigi tasks this is achieved by having all LocalTargets look like:

    .. code-block:: python

        LocalTarget(os.path.join(self.base_dir, VERSION, PIPELINE, folder, file))
        LocalTarget(os.path.join(self.scratch_dir, VERSION, PIPELINE, folder, file))


The VERSION and PIPELINE variables are defined at the top of the file:

    .. code-block:: python

        FILE_HASH = utils.file_hash(__file__)
        VERSION = fieldpathogenomics.__version__.rsplit('.', 1)[0]
        PIPELINE = os.path.basename(__file__).split('.')[0]


Maintaining structured outputs like this is crucial as luigi uses the existence of LocalTargets to determine whether or not a Task is complete.



