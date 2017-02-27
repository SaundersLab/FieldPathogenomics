import pickle
import sys
import os


def main():
    os.chdir(sys.argv[2])
    with open(sys.argv[1], 'rb') as jobfile:
        job = pickle.load(jobfile)

    os.chdir(sys.argv[3])
    if isinstance(job, dict):
        job_cls = type('job', (), job)
        job_cls().work()
    else:
        job.work()


if __name__ == '__main__':
    main()
