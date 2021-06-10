import os
import glob

def _findNumber(idx, fn):
    while idx < len(fn) and fn[idx] in "0123456789":
        idx += 1

    if idx < len(fn) and fn[idx] == ".": 
        # some other numbers must follow, to make sure it's not the start of file type suffix
        idx += 1
        digitCount = 0
        while idx < len(fn) and fn[idx] in "0123456789":
            idx += 1
            digitCount += 1

        if digitCount == 0: # if no digits follow, the '.' should not be taken into account
            idx -= 1
    
    return idx
    
def _matches(fn, fmt):

    idx = fmt.find("zZ")
    if idx < 0: # No wildcard present, look for 'z'
        idx = fmt.find("z")
        if idx < 0:
            raise Exception("Can't seem to find any redshift: neither 'zZ' nor 'z' is present")
    
        numStart = idx+1
        numEnd = _findNumber(numStart, fn)
    else:
        fmtSuffix = fmt[idx+2:]
        if fn[:idx+1] != fmt[:idx+1]:
            raise Exception("Filename doesn't match format")

        numStart = idx+1
        numEnd = _findNumber(numStart, fn)

        if fn[numEnd:] != fmtSuffix:
            return None, None
    
    return fn, float(fn[numStart:numEnd])

# using 'zZ' as indicator for redshift
def getFileNamesAndRedshifts(fileFmt, subDir = "."):
    curDir = os.getcwd()
    results = []
    try:
        os.chdir(subDir)
        for i in glob.glob(fileFmt.replace("zZ", "z*")):
            fn, z = _matches(i, fileFmt)
            if not fn:
                continue

            results.append([os.path.join(subDir, fn), z])
    finally:
        os.chdir(curDir)
    return sorted(results)
