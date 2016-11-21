import luigi
import os
class CheckTargetNonEmpty(object): 

    def complete(self):
        print("OVERRIDEN")
        outputs = luigi.task.flatten(self.output())
        return super().complete() and all(map(lambda output: os.path.getsize(output.path) > 0, outputs))


