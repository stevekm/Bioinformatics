#!/usr/bin/env python
"""
Script to get the target_capture_kit_target_region for a bam file from TCGA GDC

Example URL:
https://api.gdc.cancer.gov/files/d3a1107a-ec80-4c99-9135-f11aa0a9202d?pretty=true&expand=analysis.metadata.read_groups
https://api.gdc.cancer.gov/files/041d9457-2ffe-4d0f-8a65-1bb95e92c570?pretty=true&expand=analysis.metadata.read_groups
"""
import sys
import requests
import json
UUID = sys.argv[1]
# UUID = 'd3a1107a-ec80-4c99-9135-f11aa0a9202d'
# UUID = '041d9457-2ffe-4d0f-8a65-1bb95e92c570' # <- multiple dead URLs

output_file = sys.argv[2]

r = requests.get('https://api.gdc.cancer.gov/files/{0}?pretty=true&expand=analysis.metadata.read_groups'.format(UUID))

region_URLs = []

if r.json().get('data'):
    for read_group in r.json()['data']['analysis']['metadata']['read_groups']:
        if read_group['target_capture_kit_target_region'] is not None:
            for url in read_group['target_capture_kit_target_region'].split('|'):
                region_URLs.append(url)

d = {
'tumor_bam_uuid': UUID,
'target_capture_kit_target_regions': list(set(region_URLs))
}

with open(output_file, 'w') as f:
    json.dump(d, f, indent = 4)
