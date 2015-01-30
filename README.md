Unique-Sam
==========

**Unique-Sam** is a simple command line tool to remove the duplicated alignments in the [SAM](https://github.com/samtools/hts-specs) file. If the MAPQ field of the alignment is available, *unique-sam* will keep one and only one alignment with the highest score. Otherwise, *unique-sam* will calculate a score according to the alignment's MD or CIGAR field and use the calculated value to remove the duplicated alignments.

Install
=====
- Install with the source code, in the source folder:
```bash
python setup.py install
```

- If you have [**pip**](https://pip.pypa.io/en/latest/index.html) installed, you can simply run 
```bash
pip install unique-sam
``` 
After installation you can access **unique-sam** from your command line.

Usage
=====
**unique-sam** need a SAM format file to run properly. In your command line environment:
```bash
unique-sam input.sam -o output.sam
```
If no `-o output.sam` option is specified, the output will be written to `unique.inputfilename.sam` by default. 

For more about **unique-sam** options run:
```bash
unique-sam --help
```  
Copyright
========
Copyright (c) 2015 [dlmeduLi@163.com](mailto:dlmeduLi@163.com)


