from luigi.contrib import sqla
from luigi.task import flatten
import luigi
import os
import sqlalchemy
import datetime
import fieldpathogenomics.utils as utils


class CommitToTable(sqla.CopyToTable):
    columns = [(["path", sqlalchemy.String(4096)], {}),
               (["checksum", sqlalchemy.INTEGER], {}),
               (["datetime", sqlalchemy.DateTime], {}),
               (["task_family", sqlalchemy.String(100)], {}),
               (["git_commit", sqlalchemy.String(40)], {}),
               (["pipeline_hash", sqlalchemy.String(40)], {})]

    connection_string = "mysql+pymysql://tgac:tgac_bioinf@tgac-db1.hpccluster/buntingd_fieldpathogenomics"
    table = "FileTable"
    _rows = luigi.ListParameter()

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def rows(self):
        return self._rows

    def copy(self, conn, ins_rows, table_bound):
        """
        Override copy to first delete any stale paths
        """
        # First delete any existing entry for the same path
        for row in ins_rows:
            rm = table_bound.delete().where(table_bound.c.path == row['_path'])
            conn.execute(rm)
        super().copy(conn, ins_rows, table_bound)

    def update_id(self):
        return hash(str(self._rows))


class CommittedTarget(luigi.LocalTarget):
    '''LocalTarget that is checksummed and git commit hash stored
       in the file database on creation.
       Use this for any permenant output files, this allows us to verify their
       integrity later and also related the output file to a particular version of the
       code that created it.'''
    def checksum(self):
        return utils.checksum(self.path)

    def commit(self, task_family, pipeline_hash):
        row = (os.path.abspath(self.path),
               self.checksum(),
               str(datetime.datetime.now()),
               task_family,
               utils.current_commit_hash(os.path.split(__file__)[0]),
               pipeline_hash)
        CommitToTable([row]).run()


class CommittedTask():
    '''Any task that creates CommittedTargets need to subclass this mixin.
       It overrides on_sucess to commit the taget checksum to SQL'''
    def on_success(self):
        pipeline_hash = utils.hash_pipeline(self)
        for o in flatten(self.output()):
            if isinstance(o, CommittedTarget):
                o.commit(task_family=self.task_family,
                         pipeline_hash=pipeline_hash)
