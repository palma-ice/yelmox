{
    "defaults" :
    {
        "jobname"       : "Yelmo",
        "email"         : "USER@awi.de",
        "group"         : "envi.envi", 
        "omp"           : 0,
        "wall"          : "48:00:00", 
        "qos"           : "48h",
        "partition"     : "smp",
        "job_template"  : "config/albedo_submit_slurm"
    },

    "exe_aliases" :
        {   "yelmox" : "libyelmox/bin/yelmox.x",
            "iso"    : "libyelmox/bin/yelmox_iso.x",
            "hyst"   : "libyelmox/bin/yelmox_hyst.x",
            "rembo"  : "libyelmox/bin/yelmox_rembo.x",
            "ismip6" : "libyelmox/bin/yelmox_ismip6.x"
        },

    "grp_aliases" : {},

    "par_paths" : 
        {
            "rembo"  : "par/rembo_Greenland.nml"
        },

    "files" : ["git_yelmo.txt"], 

    "dir-special" : {},

    "links" : 
        ["input","ice_data","maps"],

    "job_queues" :
        {   "12h" :
            {   "wall" : "12:00:00"  },
            "48h" :
            {   "wall" : "48:00:00"  },
            "30min" :
            {   "wall" : "00:30:00" },
            "1wk" :
            {   "wall" : "168:00:00" }
        }
}
