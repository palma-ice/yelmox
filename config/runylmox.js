{
    "defaults" :
    {
        "jobname"       : "Yelmo",
        "email"         : "USER@pik-potsdam.de",
        "group"         : "anthroia", 
        "omp"           : 0,
        "wall"          : "24:00:00", 
        "qos"           : "priority",
        "partition"     : "haswell",
        "job_template"  : "config/pik_submit_slurm"
    },

    "exe_aliases" :
        {   "yelmox" : "libyelmox/bin/yelmox.x",
            "iso"    : "libyelmox/bin/yelmox_iso.x",
            "hyst"   : "libyelmox/bin/yelmox_hyst.x",
            "rembo"  : "libyelmox/bin/yelmox_rembo.x",
            "ismip6" : "libyelmox/bin/yelmox_ismip6.x",
            "rtip"   : "libyelmox/bin/yelmox_rtip.x",
        },

    "grp_aliases" : {},

    "par_paths" : 
        {
            "rembo"  : "par/rembo_Greenland.nml"
        },

    "files" : ["git_yelmo.txt"], 

    "dir-special" : {},

    "links" : 
        ["input","ice_data","isostasy_data","maps"],

    "job_queues" :
        {   "priority" :
            {   "wall" : "24:00:00"  },
            "short" :
            {   "wall" : "24:00:00"  },
            "medium" :
            {   "wall" : "168:00:00" },
            "long" :
            {   "wall" : "720:00:00" }
        }
}
