import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import neospy

jd_start = neospy.wise.MISSION_PHASES['Cryo'].jd_start
jd_end = neospy.wise.MISSION_PHASES['Reactivation_2023'].jd_end

# These objects are explicitly ignored
ignore = {'483P',       # This currently has record problems between JPL and MPC
          'A/2018 V3',  # TNO - Horizons doesnt classify it's orbit as a comet, MPC does
          'A/2020 A1',  # TNO - Horizons doesnt classify it's orbit as a comet, MPC does
          'C/2010 E3',  # Single epoch of observation by WISE only, bad orbit fit
          'C/2019 L2',  # Single epoch of observation by WISE only, bad orbit fit
         }

def download_comet_orbit_data(filename = "jpl_comet.csv"):
    """
    ** THIS DOESNT NEED TO BE RUN, THE FILE PRODUCED BY THIS IS ALREADY PRESENT **

    This function is left here as a record of how that file was made.
    
    Download comet orbit metadata from JPL Horizons.

    This constructs the jpl_comets.csv file which contains the information about comet orbits.

    These comet orbits are made up of 4 elements:
        1) Record ID, a unique record ID used by JPL internals to keep track of unique orbits.
        2) Epoch year, the year for the epoch of fit of the orbit.
        3) Desig, The designation for the object.
        4) Name, the full name of the object
    """
    # Downloading the set of all comets from JPL
    import requests
    res = requests.get("https://ssd.jpl.nasa.gov/api/horizons.api?command='COM%3BNOFRAG%3B'")
    
    rows = []
    for row in res.json()['result'].split("\n")[10:-4]:
        row = [row[:16].strip(), row[16:22].strip(), row[22:38].strip(), row[38:].lstrip()]
        rows.append(", ".join([f"{r:>10s}" for r in row]).strip()+"\n")
    with open(filename, "w") as f:
        f.write("# "+", ".join([f"{r:>10s}" for r in ['record','epoch', 'desig', 'name          ']]).strip()+"\n")
        f.writelines(rows)
        
def load_comet_orbit_metadata(filename = "jpl_comet.csv", keep_after=1980):
    """
    Load the comet orbit metadata, trimming the epoch by the provided year limit.

    See the above function for definitions
    """
    jpl_comets = pd.read_csv(filename, names=['record', 'epoch', 'desig', 'name'], header=1)
    jpl_comets = jpl_comets[jpl_comets.epoch > keep_after]
    jpl_comets.desig = [s.strip() for s in jpl_comets.desig]
    return jpl_comets

def load_mpc_observations():
    """
    Load observations from the minor planet center for comets.
    This trims the observation list to only include the folllowing:
        1) C51 (WISE) observations.
        2) Remove 'A/' designations.
        3) Remove objects in the ignore list at the top of this file.
        4) Remove observations after jd_end specified at the top of this file.
    """
    coms = neospy.cache.cached_file_download("https://minorplanetcenter.net/iau/ECS/MPCAT-OBS/CmtObs.txt.gz")
    import gzip
    with gzip.open(coms, 'r') as f:
        lines = [x.decode() for x in f.readlines()]
    lines = [l for l in lines if "C51" in l[-9:]]
    comet_obs = neospy.mpc.MPCObservation.from_lines(lines)
    comet_obs = [c for c in comet_obs if c.obs_code =="C51"]
    comet_obs = [c for c in comet_obs if c.jd <= jd_end]
    for obs in comet_obs:
        if len(obs.name) > 8:    
            obs.name = neospy.mpc.unpack_comet_designation(obs.name[:5])
        else:
            obs.name = neospy.mpc.unpack_comet_designation(obs.name)

        # Rename some objects which MPC names differently from JPL
        if obs.name == "P/2010 CG6":
            obs.name = "P/2010 J3"
        if obs.name == "2I":
            obs.name = "C/2019 Q4"
        if obs.name == "P/2010 PB57":
            obs.name = "P/2010 WK"
    
    comet_obs = [c for c in comet_obs if "A/" not in c.name]
    comet_obs = [c for c in comet_obs if c.name not in ignore]
    
    return comet_obs

def filter_names_using_horizons(desigs, pos_unc=0.01, min_dist=11.5):
    """
    Given a list of designations, lookup the name which horizons prefers for this object.

    This rejects objects with the following properties:
        1) If the covariance of the orbit fit is not specified (typically very bad orbit fits).
        2) If positional uncertainty is greater than `pos_unc` specified.
        3) Reject if the objects minimum distance to the sun during the WISE mission is never
           less than the specified `min_dist`

    This returns a dictionary from the JPL preferred name to the original names provided.
    """

    found_names = {}
    for idx, name in enumerate(desigs):
        if "D/" in name:
            continue
        original_name = name
        try:
            obj = neospy.HorizonsProperties.fetch(name, exact_name=True)
        except:
            prefix, name = name.split("/", maxsplit=1)
            obj = neospy.HorizonsProperties.fetch(name, exact_name=True)
        if obj.covariance is None:
            continue
        pos_cov = np.linalg.norm(np.diag(obj.covariance.cov_matrix)[:3])
        if pos_cov > pos_unc:
            continue
                
        state = obj.state
        
        d_min = 1000
        for jd in np.arange(jd_start, jd_end, 40):
            if neospy.wise.mission_phase_from_jd(jd) is None:
                continue
            d = neospy.propagate_two_body([state], jd)[0].pos
            d = (neospy.SpiceKernels.state("Earth", jd).pos - d).r
            if d < d_min:
                d_min = d
        if d_min > min_dist:
            continue
                
        found_names[name] = original_name
    return found_names

def load_fovs():
    """
    Load a full set of the field of views from WISE.
    
    **  This takes more than 50GB of ram.  **
    
    """
    fovs = []
    for phase in neospy.wise.MISSION_PHASES.values():
        f = list(neospy.wise.fetch_WISE_fovs(phase))
        fovs.extend(f)
    return fovs

def calculate_vis():
    """
    ** THIS DOESNT NEED TO BE RUN, THE FILE PRODUCED BY THIS IS ALREADY PRESENT **

    This function is left here as a record of how that file was made.
    
    Given all SPK objects available in the kernel folder, calculate all visible objects
    which are comet designations.

    This takes several hours on a desktop and at least 50 GB of ram.

    Results are saved to a master `wise_all.bin` file.
    """
    neospy.SpiceKernels.cache_kernel_reload()
    fovs = load_fovs()
    spkids = [s[1] for s in neospy.SpiceKernels.loaded_objects() if s[1] > 1000 and s[1] < 20000001]
    vis = [v[0] for v in neospy.propagation.spice_visible(spkids, fovs)]
    neospy.SimultaneousStates.save_list(vis, f"wise_all.bin")
    return vis

def load_vis():
    """
    Load the pre
    """
    import tarfile
    from tempfile import TemporaryDirectory
    
    with TemporaryDirectory() as t:
        with tarfile.open("wise_all.tar.xz") as f:
            f.extractall(path=t)
            vis = neospy.SimultaneousStates.load_list(t + "/wise_all.bin")
    return vis


def download_spks(keep):
    """    
    ** THIS DOESNT NEED TO BE RUN, FINAL DATA IS ALREADY PROVIDED **

    This function is left here as a record of how that file was made.
    
    Download all spice kernels for objects in the "keep" list.
    """
    final_jd_start = neospy.Time.from_ymd(2000, 1, 1).jd
    final_jd_end = neospy.Time.from_ymd(2035, 1, 1).jd

    
    found_names = filter_names_using_horizons(keep)
    found_names_rev = {v: k for k, v in found_names.items()}
    
    for idx, k in enumerate(sorted(keep)[:]):
        # This comet needs to be manually requested from horizons, since it is an odd-ball
        if k == "C/2013 A1":
            continue

        subset = jpl_comets[jpl_comets.desig == k]
        assert len(subset) > 0
        epochs = np.array(subset.epoch)
        
        mid_points = (epochs[0:-1] + epochs[1:]) / 2
        mid_points = [neospy.Time.from_ymd(m, 6, 1).jd for m in mid_points]
        
        
        segments = []
        jd_end = final_jd_end
        for step in mid_points[::-1]:
            if step < final_jd_start:
                break
            segments.append([step, jd_end])
            jd_end = step
            
        segments.append([final_jd_start, jd_end])
        segments.reverse()
        segments
    
        for seg in segments:
            print(idx, "/", len(keep), "   - ", k, " Dates: ", neospy.Time(seg[0]).ymd[0], "-", neospy.Time(seg[1]).ymd[0])
            for _ in range(5):
                try:
                    neospy.SpiceKernels.cached_kernel_horizons_download(found_names_rev[k], seg[0], seg[1], exact_name=True)
                    break
                except Exception as e:
                    print(e)
                    
def make_keep_list():
    """
    ** THIS DOESNT NEED TO BE RUN, THE FILE PRODUCED BY THIS IS ALREADY PRESENT **

    This function is left here as a record of how that file was made.
    
    Create the list of object names which are kept under the filter criterion.
    """
    found_names = filter_names_using_horizons(desigs)
    found_names_rev = {v: k for k, v in found_names.items()}
    
    with open("keepers.txt", "w") as f:
        f.write("\n".join(sorted(list(found_names_rev.keys()))))
    
    keep = sorted(list(found_names_rev.keys()))
    return keep

def load_keep_list():
    """
    Load the list of object names which are kept under the filter criterion.
    """
    with open("keepers.txt", 'r') as f:
        return [f.strip() for f in f.readlines()]
