{
    "defaults" :
    {
        "jobname"       : "Yelmo",
        "email"         : "USER@ucm.es",
        "group"         : "", 
        "omp"           : 0,
        "wall"          : 24, 
        "qos"           : "normal",
        "partition"     : "short",
        "job_template"  : "config/brigit_submit_slurm"
    },

    "exe_aliases" :
        {   "yelmox" : "libyelmox/bin/yelmox.x",
            "iso"    : "libyelmox/bin/yelmox_iso.x",
            "hyst"   : "libyelmox/bin/yelmox_hyst.x",
            "rembo"  : "libyelmox/bin/yelmox_rembo.x"
        },

    "grp_aliases" : {},

    "par_paths" : {},

    "files" : ["par/rembo_Greenland.nml","git_yelmo.txt"], 

    "dir-special" : {},

    "links" : 
        ["input","ice_data","maps"],

    "job_queues" :
        {  "normal" :
            {   "wall" : 10000 }
        }
}
