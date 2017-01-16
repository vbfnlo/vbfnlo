#!/usr/bin/env python3
import re


def replacemachine(fileName, sourceText, replaceText, keep=True):
    with open(fileName, "r") as file:
        text = file.read()
    with open(fileName, "w") as file:
        if keep:
            file.write(re.sub(r'^[\s!]*'+sourceText+r'\s*=',
                              sourceText+' = '+replaceText+' ! ',
                              text,
                              flags=re.MULTILINE))
        else:
            file.write(re.sub(r'^[\s!]*'+sourceText+r'\s*=.*',
                              sourceText+' = '+replaceText+' \n',
                              text,
                              flags=re.MULTILINE))


def readconf(filename):
    COMMENT_CHAR = '!'
    OPTION_CHAR = '='
    options = {}
    with open(filename) as f:
        for line in f:
            # First, remove comments:
            if COMMENT_CHAR in line:
                # split on comment char, keep only the part before
                line, comment = line.split(COMMENT_CHAR, 1)
            # Second, find lines with an option=value:
            if OPTION_CHAR in line:
                # split on option char:
                option, value = line.split(OPTION_CHAR, 1)
                # strip spaces:
                option = option.strip()
                value = value.strip()
                # store in dictionary:
                options[option] = value
    return options
