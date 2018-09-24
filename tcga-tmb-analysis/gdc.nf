params.manifest = "gdc_manifest.2018-09-13.txt"
params.outputDir = "."
params.server = 'https://api.gdc.cancer.gov/'

process gdc_download {
    echo true
    executor "local"
    storeDir "${params.outputDir}"

    input:
    file(manifest) from Channel.fromPath("${params.manifest}")

    output:
    file('gdc')

    script:
    """
    which gdc-client
    gdc-client --version

    mkdir gdc
    gdc-client download -m "${manifest}" -s "${params.server}" -d gdc -n 4
    for item in \$(find gdc/ -type f -name "*.gz"); do
    echo ">>> Unzipping \${item}"
    gunzip "\${item}"
    done

    # wget https://genome.nyumc.org/results/external/NYU/snuderllab/data/TCGA/gdc.tar.gz
    # tar -xzf gdc.tar.gz
    # rm -f gdc.tar.gz
    """
}
