import pickle
import sys
import os


def main():
    os.chdir(sys.argv[2])
    with open(sys.argv[1], 'rb') as jobfile:
        job = pickle.load(jobfile)
    os.chdir(sys.argv[3])
    job.work()


if __name__ == '__main__':
    main()
