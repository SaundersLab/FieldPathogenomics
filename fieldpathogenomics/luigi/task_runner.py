import dill
import sys


def main():
    with open(sys.argv[1], 'rb') as jobfile:
        job = dill.load(jobfile)
    job.work()


if __name__ == '__main__':
    main()
