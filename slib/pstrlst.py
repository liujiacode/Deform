#!/usr/bin/env python3
# Paraphrasing string and list.
#
# Version 1.0.0
# 2019.01.13
#
# Author: Liu Jia
#


import re


# Delete head and tail spaces in string (stri).
def del_stri_ht_space(stri, del_chars=(" ", "\t", "\n", "\r")):
    if not isinstance(stri, str):
        raise ValueError("expect str for stri.")
    strii = stri
    while strii and strii[0] in del_chars:
        strii = strii[1:]
    while strii and strii[-1] in del_chars:
        strii = strii[:-1]
    return strii

# Delete all the spaces in list (lsti).
def del_lsti_al_space(lsti, del_chars=(" ", "\t", "\n", "\r")):
    if not isinstance(lsti, list):
        raise ValueError("expect list for lsti.")
    lstii = lsti
    for i in range(len(lstii)):
        if isinstance(lstii[i], str) and del_stri_ht_space(lstii[i], del_chars=del_chars) == "":
            lstii[i] = ""
    while "" in lstii:
        lstii.remove("")
    return lstii

def del_comment(stri, signal=("#"), del_chars=(" ", "\t", "\n", "\r")):
    strii = del_stri_ht_space(stri, del_chars=del_chars)
    line = re.split("|".join(signal), strii)[0]
    line = del_stri_ht_space(line, del_chars=del_chars)
    return line


if __name__ == "__main__":
    # Test cases.
    assert del_comment("               ") == ""
    assert del_comment("# abc = 123 456") == ""
    assert del_comment("abc # = 123 456") == "abc"
    assert del_comment("abc = # 123 456") == "abc ="
    assert del_comment("abc = 123 # 456") == "abc = 123"
    assert del_comment("abc = 123 456 #") == "abc = 123 456"
    assert del_comment("  #              ") == ""
    assert del_comment("  # abc = 123 456") == ""
    assert del_comment("  abc # = 123 456") == "abc"
    assert del_comment("  abc = # 123 456") == "abc ="
    assert del_comment("  abc = 123 # 456") == "abc = 123"
    assert del_comment("  abc = 123 456 #") == "abc = 123 456"
    assert del_comment("                 #") == ""
    assert del_comment("# abc = 123 456   ") == ""
    assert del_comment("abc # = 123 456   ") == "abc"
    assert del_comment("abc = # 123 456   ") == "abc ="
    assert del_comment("abc = 123 # 456   ") == "abc = 123"
    assert del_comment("abc = 123 456 #   ") == "abc = 123 456"
