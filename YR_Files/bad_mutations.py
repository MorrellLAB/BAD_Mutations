#!/usr/bin/env python
import argparse

def parse_args():
    parser = argparse.ArgumentParser(description='Process some integers.')
    parser.add_argument('--user', required=True, help='user name for jgi.doe.gov')
    parser.add_argument('--password', required=True, help='password for jgi.doe.gov')
    args = parser.parse_args()
    return args

def main():
    print "main"

main()
