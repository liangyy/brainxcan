import os
import workflow
run_snmk = os.path.join(workflow.basedir, "brainxcan/snmk/run.snmk")
include: run_snmk
        