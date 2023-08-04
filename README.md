# Pedify Toolkit: Algorithm And Commend line program For Crop Pedigree Identify And Reconstruct.


## Dependence

This project was based on python3 environment. If your operating system does not have Python3 environment installed, Please install the [Python3](https://www.python.org/downloads/) first.


## Quickstart

### Step1: clone this code


> git clone https://github.com/plantdna/pedify.git

#### OR

> [download from here](https://github.com/plantdna/pedify/archive/refs/heads/main.zip)


### Step2: install dependent packages and activate virtual environment

``` bash
cd pedify toolkit
python3 -m venv venv
pip3 install -r requirements.txt
source /venv/bin/activate
```

### Step3: run the pendify.py

```bash
python3 pedify.py -c "config.json"
```

## config.json
> This is the parameter profile of the program, the implications and types of related parameters are shown:

- mode
    - BD : build database
    - PI : pedigree identify
    - PR : pedigree reconstruct
    - DPR : derived lines pedigree reconstruct

|key|type|meaning|required|
|---|----|-------|--------|
|mode|String|analysis mode include ["BD", "PI", "PR", "DPR"]|true|
|locus_file|String|file path of locus information|true|
|HTP|Boolean|whether to use HTP include [true, false]|false|
|output|String|output file path|true|
|miss|String|missing data default "---"|false|
|cores|Int|number of CPU cores default 1|false|
|sep|String|file separator default "\t"|false|
|param|Dict||true|

### param list if mode is BD
|key|type|meaning|required|
|---|----|-------|--------|
|genotypes|String|folder path of sample genotype files|true|
|maf_limit|Float|default 0.05|false|
|pic_limit|Float|default 0.02|false|

### param list if mode is PI
|key|type|meaning|required|
|---|----|-------|--------|
|target|String|genotype file path of target|true|
|ancestors|List|genotype files path of ancestors|true|
|dataset|String|file path of dataset|true|
|algorithm|String|the algorithm of predigree identification include ["WPI", "LPI"]|true|

### param list if mode is PR
|key|type|meaning|required|
|---|----|-------|--------|
|target|String|genotype file path of target|true|
|ancestors|Dict|genotype files path of ancestors suct as {"U8112": "./files/u8112.txt"}|true|
|colors|Dict|set ancestors color suct as {"U8112": "#fa8c16"}|true|

### param list if mode is DPR
|key|type|meaning|required|
|---|----|-------|--------|
|ancestor|String|genotype file path of ancestor|true|
|ancestorName|String|name of ancestor|true|
|derivedlines|Dict|genotype files path of derivedlines suct as {"U8112": "./files/u8112.txt"}|true|
|colors|Dict|set ancestor and derived color suct as {"ancestor": "#16C2C1", "derived": "#1777FF"}|true|


## File Examples

### sample genotype file

> missing data use "---" by default; use "-" represent missing for InDel locus

|locus|gt|
|---|----|
|locus1|A/A|
|locus2|A/T|
|locus3|---|
|locus4|-/A|

### locus information file

|locus|chr_id|position|HTP|
|-----|------|--------|---|
|locus1|1|3497|HTP0001|
|...|...|...|....|
|locusN|10|149576018|HTP6163|



## config.json Examples

### Build the dataset

```json
{
    "mode": "BD",
    "param": {
        "genotypes": "folder path",
        "maf_limit": 0.05,
        "pic_limit": 0.02
    },
    "output": "folder path",
    "locus_file": "./locus_info.txt",
    "sep": "\t"
}
```

### Pedigree Identify

```json
{
    "mode": "PI",
    "param": {
        "target": "file path",
        "ancestors": [
            "file path"
        ],
        "dataset": "file path",
        "algorithm": "WPI"
    },
    "output": "folder path",
    "locus_file": "./locus_info.txt",
    "sep": "\t",
    "HTP": true,
    "miss": "---",
    "cores": 5
}
```

### Pedigree Reconstruct

```json
{
    "mode": "PR",
    "param": {
        "target": "file path",
        "ancestors": {
            "key1": "file path",
            "key2": "file path",
        },
        "colors": {
            "key1": "#16C2C1",
            "key2": "#1777FF",
        }
    },
    "output": "folder path",
    "locus_file": "./locus_info.txt",
    "sep": "\t",
    "HTP": true,
    "window_size": 10,
    "species":"maize",
    "miss": "---"
}
```

### Derived lines Pedigree Reconstruct

```json
{
    "mode": "DPR",
    "param": {
        "ancestor": "file path",
        "ancestorName": "name",
        "derivedlines": {
            "key1": "file path",
            "key2": "file path",
        },
        "colors": {
            "ancestor": "#16C2C1",
            "derived": "#1777FF"
        }
    },
    "output": "folder path",
    "locus_file": "./dataset/locus_info.txt",
    "sep": "\t",
    "HTP": true,
    "window_size": 10,
    "miss": "---"
}
```
