<p align="center">
  <img src="https://raw.githubusercontent.com/LiuzLab/AI_MARRVEL/main/docs/images/logo.v1.png" alt="AI-MARRVEL logo"/>
</p>

# AI-MARRVEL

<p align="center">
  <a href="#" style="text-decoration:none">
    <img src="https://img.shields.io/badge/AI_MARRVEL-v1.0.1-blue.svg" alt="AI_MARRVEL version badge"/>
  </a>
  <a href='https://ai-marrvel.readthedocs.io/en/latest/?badge=latest'>
      <img src='https://readthedocs.org/projects/ai-marrvel/badge/?version=latest' alt='Documentation Status' />
  </a>
</p>

**AI-MARRVEL (AIM)** is an AI system for rare genetic disease diagnosis.

It takes as input patient VCF and phenotype (formatted with HPO) to predict the causal variant(s).  
In making prediction, it takes variant annotation from [MARRVEL](https://marrvel.org/) database and more, 
and generates **prediction score** + **confidence score** as output.

You can use AI-MARRVEL from our [website](https://ai.marrvel.org/) or follow the [documentation](https://ai-marrvel.readthedocs.io/en/main/) to run locally.

:new: Our paper is now published in [NEJM AI](https://ai.nejm.org/doi/full/10.1056/AIoa2300009)!

---

## Quick Start

### Install Required Data Dependencies

AIM utilizes various databases for variant annotation, all of which have been compiled and are available for download. We use AWS S3 for data access, and the data can be downloaded by following these steps:

1. **Install the AWS CLI**: Follow the instructions provided in the [AWS CLI Installation Guide](https://docs.aws.amazon.com/cli/latest/userguide/getting-started-install.html).
2. **Navigate to Your Desired Directory**: Change to the directory where you want your data dependencies downloaded. For example, in Ubuntu, use:
3. Use the following command to sync the S3 bucket to your local directory:

```bash
aws s3 sync s3://aim-data-dependencies-2.0-public . --no-sign-request
```

### Get the software

AIM is released as a Nextflow pipeline for easy distribution. To get it:

```bash
git clone https://github.com/LiuzLab/AI_MARRVEL
cd AI_MARRVEL
nextflow run main.nf --version
```

### Run with your sample

```bash
nextflow run main.nf  --ref_dir <PATH_TO_REFERENCE_DIRECTORY>
                  --input_vcf <PATH_TO_INPUT_VCF_FILE>
                  --input_hpo <PATH_TO_INPUT_HPO_FILE>
                  --outdir <PATH_TO_OUTPUT_DIRECTORY>
                  --bed_filter <PATH_TO_BED_FILE> # Optional
                  --run_id [Sample Id] # Optional, default: 1
                  --ref_ver [Reference genome: hg19/hg38] # Optional, default: hg19
                  --exome_filter # Optional
```

Alternatively, the pipeline can be executed with a parameter file (yaml)

```bash
nextflow run main.nf -params-file params.yaml
```

NOTE: You need to create `params.yaml` by copying [params.yaml.example](params.yaml.example) file and follow the instruction.

For more information on usage and parameters which are open for modification, please use `--help` option as shown below.

```

nextflow run main.nf --help
```

## License

AI-MARRVEL is licensed under GPL-3.0. You are welcomed to use it for research purpose.  
For business purpose, please contact us for licensing.

## Disclaimer

- Some of the data and software included in the distribution may be subject to third-party constraints. Users of the data and software are solely responsible for establishing the nature of and complying with any such restrictions.
- AI-MARRVEL provides this data and software in good faith, but make no warranty, express or implied, nor assume any legal liability or responsibility for any purpose for which they are used.

## Citing AI-MARRVEL

```bibtex
@article{doi:10.1056/AIoa2300009,
author = {Dongxue Mao  and Chaozhong Liu  and Linhua Wang  and Rami AI-Ouran  and Cole Deisseroth  and Sasidhar Pasupuleti  and Seon Young Kim  and Lucian Li  and Jill A. Rosenfeld  and Linyan Meng  and Lindsay C. Burrage  and Michael F. Wangler  and Shinya Yamamoto  and Michael Santana  and Victor Perez  and Priyank Shukla  and Christine M. Eng  and Brendan Lee  and Bo Yuan  and Fan Xia  and Hugo J. Bellen  and Pengfei Liu  and Zhandong Liu },
title = {AI-MARRVEL — A Knowledge-Driven AI System for Diagnosing Mendelian Disorders},
journal = {NEJM AI},
volume = {1},
number = {5},
pages = {AIoa2300009},
year = {2024},
doi = {10.1056/AIoa2300009},

URL = {https://ai.nejm.org/doi/abs/10.1056/AIoa2300009},
eprint = {https://ai.nejm.org/doi/pdf/10.1056/AIoa2300009}
,
    abstract = { AI-MARRVEL is an AI system for genetic diagnosis that improves diagnostic accuracy, surpassing state-of-the-art benchmarked methods. }
}
```
