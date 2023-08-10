{
  "format_mito_ref": {
    "name": "{rule}_mitogenome",
    "qos": "test",
    "time": "00:05:00",
    "queue": "genD",
    "threads": 4
  },
  "align_ont_to_mitoref": {
    "name": "{rule}_mitogenome",
    "qos": "normal",
    "time": "6:00:00",
    "queue": "genD",
    "mem": "50G",
    "threads": 16
  },  
  "align_illumina": {
    "name": "{rule}_mitogenome",
    "qos": "normal",
    "time": "6:00:00",
    "queue": "genD",
    "mem": "50G",
    "threads": 16
  },
  "flye_assembly_meta": {
    "name": "{rule}_mitogenome",
    "qos": "normal",
    "time": "12:00:00",
    "queue": "genD",
    "mem": "100G",
    "threads": 16
  },
  "select_illumina_from_circular": {
    "name": "{rule}_mitogenome",
    "qos": "normal",
    "time": "1:00:00",
    "queue": "genD",
    "mem": "20G",
    "threads": 4
  },
  "nextpolish_circular_contig": {
    "name": "{rule}_mitogenome",
    "qos": "short",
    "time": "02:00:00",
    "queue": "genD",
    "mem": "30G",
    "threads": 16
  },
  "orient_polished_contig": {
    "name": "{rule}_mitogenome",
    "qos": "short",
    "time": "02:00:00",
    "queue": "genD",
    "mem": "12G",
    "threads": 2
  },
  "orient_reference_MT": {
    "name": "{rule}_mitogenome",
    "qos": "short",
    "time": "02:00:00",
    "queue": "genD",
    "mem": "12G",
    "threads": 2
  }
  
}
