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


I find the easiest way to access the database is to use a Jupyter notebook and pandas to query. For example to plot the percentage of mapped reads per library:

.. code-block:: python

    import pandas as pd
    import matplotlib.pyplot as plt
    import seaborn as sns
    sns.set_style('whitegrid')
    %matplotlib inline
    connection_string = "mysql+pymysql://tgac:tgac_bioinf@tgac-db1.hpccluster/buntingd_fieldpathogenomics"
    df = pd.read_sql('AlignmentStats', connection_string).apply(pd.to_numeric, args=('ignore',))
    df.plot(kind='bar', x='Library', y='mapped_reads_pc')


