import pickle
import sys


def main():
    jobfile = sys.argv[1]
    job = pickle.load(jobfile)
    job.work()


if __name__ == '__main__':
    main()
