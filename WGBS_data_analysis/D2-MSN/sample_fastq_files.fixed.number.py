#!/home/yli/anaconda3/bin/python3

import sys
import gzip

def simple_sampling(file1, file2 ):
    record_no = 0
    with gzip.open(file1) as input_file:
        with gzip.open(file2,'w') as output_file:
            for line1 in input_file:
                line2 = next(input_file)
                line3 = next(input_file)
                line4 = next(input_file)
                if (record_no % 4 == 0 )&( record_no < 4000000):
                    output_file.write(line1)
                    output_file.write(line2)
                    output_file.write(line3)
                    output_file.write(line4)
                record_no += 1
    return 0

if __name__ == '__main__':
    print("usage: python3 spl.py file1 file3 ")
    if (sys.argv[1] and sys.argv[2]):
        simple_sampling(sys.argv[1], sys.argv[2])
    else:
        raise NameError