import pickle
import sys


def main():
    with open(sys.argv[1], 'rb') as jobfile:
        job = pickle.load(jobfile)
    job.work()


if __name__ == '__main__':
    main()
