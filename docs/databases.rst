Databases
==============

MySQL server is tgac-db1.nbi.ac.uk at the moment the database is buntingd_fieldpathgenomics.

Sample Metadata Table
--------------------------
This needs to happen.
I currently have about 6 different versions of the excel files with sample information in.
I tried to unify them but it's non-trivial as there is no standardisation.
Diane/Antoine please advise!


Alignments Statistics Table
----------------------------

This table stores the STAR alignments logs and is populated by the :py:class:`Library.AlignmentStats` task.
Schema currently is 

.. code-block:: python

    columns = [
        (["Library", sqlalchemy.String(64)], {}),
        (["input_reads", sqlalchemy.INTEGER], {}),
        (["input_len", sqlalchemy.FLOAT], {}),
        (["mapped_reads", sqlalchemy.INTEGER], {}),
        (["mapped_reads_pc", sqlalchemy.String(10)], {}),
        (["mapped_len", sqlalchemy.FLOAT], {}),
        (["mismatch_pc", sqlalchemy.String(10)], {}),
        (["datetime", sqlalchemy.String(25)], {}), 
        (["genome", sqlalchemy.String(25)], {}),
        (["git_commit", sqlalchemy.String(40)], {}),
        (["pipeline_hash", sqlalchemy.String(40)], {}),
    ]


