from collections import defaultdict
import io
import sys

def get_read_info(line):
    line = line.split("_")
    return {
        'rid': line[0],
        'dir': line[1],
        'sid': int(line[2]),
        'pos': int(line[3])
            }

def main():
    # take one read
    fname = sys.argv[1]
    line_num = -1

    hist = defaultdict(int)

    with open(fname, "r") as fhandle:
        for line in fhandle:
            # skip all lines other than the ones containing sim info
            line_num += 1
            if line_num % 4 != 0:
                continue
            read_info = get_read_info(line)
            hist[read_info['sid']] += 1

    sorted_keys = sorted(hist.keys())
    vals = []
    for k in sorted_keys:
        print k, "\t", hist[k]


if __name__ == '__main__':
    main()
