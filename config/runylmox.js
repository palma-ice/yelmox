{
    "defaults" :
    {
        "jobname"       : "Yelmo",
        "email"         : "USER@pik-potsdam.de",
        "group"         : "anthroia", 
        "omp"           : 0,
        "wall"          : 24, 
        "qos"           : "priority",
        "partition"     : "haswell",
        "job_template"  : "config/pik_submit_slurm"
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
        {   "priority" :
            {   "wall" : 24  },
            "short" :
            {   "wall" : 24  },
            "medium" :
            {   "wall" : 168 },
            "long" :
            {   "wall" : 720 }
        }
}
