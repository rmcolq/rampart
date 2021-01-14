import argparse
import csv
import random

def parse_args():
    parser = argparse.ArgumentParser(description='Parse CSV of variant calls.')

    parser.add_argument("--in-file", action="store", type=str, dest="in_file")
    parser.add_argument("--out-file", action="store", type=str, dest="out_file")
    parser.add_argument("--random", action="store_true", dest="simulate")

    return parser.parse_args()

def call_mutation(key, call, simulate=False):
    if simulate:
        r = random.randint(0, 4)
        if r > 0:
            return key
        else:
            return None
    if call == 'del':
        return key
    elif call == 'X' or call == '-':
        return None
    else:
        mut = key.split(":")[-1]
        if mut.startswith(call):
            return None
        elif key.endswith(call):
            return key
        else:
            print("ERROR")
            print(key, call)
            return None

def summarize_calls(infile, outfile, simulate=False):
    calls = []
    irrelevant = ['query', 'ref_count', 'alt_count', 'other_count', 'fraction_alt']

    with open(infile) as csv_handle:
        reader = csv.DictReader(csv_handle)
        mutations = [r for r in reader.fieldnames if r not in irrelevant]
        print(mutations)

        for read in reader:
            read_name = read['query']
            read_calls = []
            for mut in mutations:
                call = call_mutation(mut, read[mut], simulate)
                if call:
                    read_calls.append(call)
            variants = "|".join(read_calls)
            calls.append({'read_name':read_name, 'variants':variants})

    with open(outfile, 'w') as out_handle:
        f = csv.DictWriter(out_handle, fieldnames=['read_name','variants'])
        f.writeheader()
        f.writerows(calls)

if __name__ == '__main__':
    args = parse_args()
    summarize_calls(args.in_file, args.out_file, args.simulate)