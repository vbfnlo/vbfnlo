#!/usr/bin/env python3
"""Functions to read and write VBFNLO-Style .dat files"""
import re

def replacemachine(filename, sourcetext, replacetext, keep=True):
    "Set value for a variable in .dat file"
    with open(filename, "r") as file:
        text = file.read()
    with open(filename, "w") as file:
        if keep:
            file.write(re.sub(r'^[\s!]*'+sourcetext+r'\s*=',
                              sourcetext+' = '+replacetext+' ! ',
                              text,
                              flags=re.MULTILINE))
        else:
            file.write(re.sub(r'^[\s!]*'+sourcetext+r'\s*=.*',
                              sourcetext+' = '+replacetext+' \n',
                              text,
                              flags=re.MULTILINE))


def readconf(filename):
    "Read config from .dat file"
    comment_char = '!'
    option_char = '='
    options = {}
    with open(filename) as file:
        for line in file:
            if comment_char in line:
                line = line.split(comment_char, 1)[0]

            if option_char in line:
                option, value = line.split(option_char, 1)
                options[option.strip()] = value.strip()
    return options
