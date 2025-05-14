import os 
import sys
import shutil
import numpy as np 
import pandas as pd
import requests
from requests.exceptions import RequestException
from zipfile import ZipFile
from tqdm import tqdm
from ase.io import read, iread



__version__ = 0.1



def available_datasets():
    """
    A dictionary containing dataset names as keys and their respective download links as values.
    """

    datasets = {
        "BVEL13k": 'https://bs3u.obs.ru-moscow-1.hc.sbercloud.ru/litraj/BVEL13k.zip',
        "nebDFT2k": 'https://bs3u.obs.ru-moscow-1.hc.sbercloud.ru/litraj/nebDFT2k.zip',
        "nebBVSE122k": 'https://bs3u.obs.ru-moscow-1.hc.sbercloud.ru/litraj/nebBVSE122k.zip',
        "MPLiTrj_subsample": 'https://bs3u.obs.ru-moscow-1.hc.sbercloud.ru/litraj/MPLiTrj_subsample.zip',
        "MPLiTrj": 'https://bs3u.obs.ru-moscow-1.hc.sbercloud.ru/litraj/MPLiTrj.zip',
        "MPLiTrj_raw": 'https://bs3u.obs.ru-moscow-1.hc.sbercloud.ru/litraj/MPLiTrj_raw.zip',
    }
    return datasets



def download_dataset(dataset_name, folder, unzip = True, remove_zip = False):
    """
    Download a specified dataset

    Parameters
    ----------

    dataset_name: str
        dataset name, can be "BVEL13k", "nebBVSE122k", "nebDFT2k", "MPLiTrj_subsample", "MPLiTrj"

    folder: str or None, None by default
        path to save a .zip file 

    unzip: boolean, True by default
        unzip donwloaded .zip
    
    remove_zip: boolean, False by default
        delete .zip file
    
    Examples
    --------

    >>> from litraj.data import download_dataset
    >>> download_dataset('nebDFT2k', '.') # download the data to the current folder

    """  

    try:
        if dataset_name not in available_datasets().keys():
            text = f'Wrong dataset name: {dataset_name}\nAllowed names: ' + ", ".join(available_datasets().keys())
            raise ValueError(text)
        
        pth = f'{folder}/{dataset_name}' if folder else f'./{dataset_name}'

        if os.path.exists(f'{pth}.zip'):
            print(f'File {pth}.zip already exists')
        else:
            _download(available_datasets()[dataset_name], filename=f'{pth}.zip')

        if unzip:
            if os.path.exists(f'{pth}'):
                print(f'Folder {pth} already exists')
            else:
                _unzip(pth, folder)
        if remove_zip:
            if os.path.exists(f'{pth}.zip'):
                os.remove(f'{pth}.zip')
                print(f"{pth}.zip has been deleted.")
                

    except KeyboardInterrupt:
        print('\nDownload interrupted by user. Cleaning up...')
        
        if os.path.exists(f'{pth}.zip'):
            os.remove(f'{pth}.zip')
        
        if unzip and os.path.exists(f'{pth}'):
            shutil.rmtree(f'{pth}')



def load_data(dataset_name, folder):
    """
    Returns a dataset

    Parameters
    ----------

    dataset_name: str
        dataset name, can be "BVEL13k", "nebBVSE122k", "nebDFT2k", "MPLiTrj_subsample", "MPLiTrj"

    folder: str
        path to the data folder

    Returns
    -------
    dataset: format depends on the dataset
        see examples

    Examples
    --------

    # nebDFT2k dataset
    >>> from litraj.data import load_data
    >>> data = load_data('nebDFT2k', 'pth/to/saved/data')
    >>> data_train = data[data._split == 'train']
    >>> for atoms in data_train.centroid:
    >>>     edge_id = atoms.info['edge_id']
    >>>     mp_id = atoms.info['material_id']
    >>>     em = atoms.info['em']
    >>>     # do stuff

    >>> from litraj.data import load_data
    >>> data = load_data('nebDFT2k', 'pth/to/saved/data')
    >>> for traj_init, traj_relaxed in zip(data.trajectory_init, data.trajectory_relaxed):
    >>>     # do stuff
    >>>     pass

    
    # MPLiTrj dataset
    >>> from litraj.data import load_data
    >>> atoms_list_train, atoms_list_val, atoms_list_test = load_data('MPLiTrj', 'pth/to/saved/data')
    >>> for atoms in atoms_list_train:
    >>>     forces = atoms.get_forces()
    >>>     energy = atoms.get_potential_energy()
    >>>     stress = atoms.get_stress()
    >>>   # do stuff

    
    # nebBVSE122k dataset
    >>> from litraj.data import load_data
    >>> atoms_list_train, atoms_list_val, atoms_list_test, index = load_data('nebBVSE122k', 'pth/to/saved/data')
    >>> for atoms_with_centroid in atoms_list_train:
    >>>     edge_id = atoms_with_centroid.info['edge_id']   # mp-id_source_target_offsetx_offsety_offsetz
    >>>     mp_id = atoms_with_centroid.info['material_id']
    >>>     em = atoms_with_centroid.info['em']
    >>>     centroid_index = np.argwhere(atoms_with_centroid.symbols =='X')
    >>>     # do stuff


    # BVEL13k dataset
    >>> from litraj.data import load_data
    >>> atoms_list_train, atoms_list_val, atoms_list_test, index = load_data('BVEL13k', 'pth/to/saved/data')
    >>> for atoms in atoms_list_train:
    >>>   mp_id =  atoms.info['material_id']
    >>>   e1d, e2d, e3d = atoms.info['E_1D'], atoms.info['E_2D'], atoms.info['E_3D']
    >>>   # do stuff

    """

    if dataset_name not in available_datasets().keys():
        text = f'Wrong dataset name: {dataset_name}\nAllowed names: ' + ", ".join(available_datasets().keys())
        raise ValueError(text)

    folder = f'{folder}/{dataset_name}'

    if dataset_name == 'nebDFT2k':
        index = pd.read_csv(f'{folder}/{dataset_name}_index.csv')
        trajectories_init, trajectories_relaxed = [], []
        for edge_id in tqdm(index.edge_id, desc = 'loading trajectories'):
            images_init = read(f'{folder}/{edge_id}_init.xyz', index =':')
            trajectories_init.append(images_init)
            images_relaxed = read(f'{folder}/{edge_id}_relaxed.xyz', index =':')
            trajectories_relaxed.append(images_relaxed)
        index['trajectory_init'] = trajectories_init
        index['trajectory_relaxed'] = trajectories_relaxed
        centroids = []
        for centroid in tqdm(iread(f'{folder}/{dataset_name}_centroids.xyz', index = ':', format = 'extxyz'),
                       desc = 'loading centroids', total = len(trajectories_init)):
            centroids.append(centroid)
        index['centroid'] = centroids
        return index
    

    if dataset_name == 'MPLiTrj_subsample':
        train, val, test = [], [], []
        for atoms in tqdm(iread(f'{folder}/{dataset_name}_train.xyz', index = ':', format = 'extxyz'),
                            desc = 'loading train', total = 94044):
            train.append(atoms)

        for atoms in tqdm(iread(f'{folder}/{dataset_name}_val.xyz', index = ':', format = 'extxyz'),
                            desc = 'loading val', total = 12249):
            val.append(atoms)

        for atoms in tqdm(iread(f'{folder}/{dataset_name}_test.xyz', index = ':', format = 'extxyz'),
                            desc = 'loading test', total = 11731):
            test.append(atoms)
        return train, val, test
    

    if dataset_name == 'MPLiTrj':
        train, val, test = [], [], []
        for atoms in tqdm(iread(f'{folder}/{dataset_name}_train.xyz', index = ':', format = 'extxyz'),
                            desc = 'loading train', total = 766254):
            train.append(atoms)

        for atoms in tqdm(iread(f'{folder}/{dataset_name}_val.xyz', index = ':', format = 'extxyz'),
                            desc = 'loading val', total = 80320):
            val.append(atoms)

        for atoms in tqdm(iread(f'{folder}/{dataset_name}_test.xyz', index = ':', format = 'extxyz'),
                            desc = 'loading test', total = 82492):
            test.append(atoms)
        return train, val, test



    if dataset_name == 'BVEL13k':
        index = pd.read_csv(f'{folder}/{dataset_name}_index.csv')
        atoms_list_train, atoms_list_val, atoms_list_test = [], [], []
        for a in tqdm(iread(f'{folder}/{dataset_name}_train.xyz', index = ':', format = 'extxyz'),
                       desc = 'loading train', total = 10159):
            atoms_list_train.append(a)


        for a in tqdm(iread(f'{folder}/{dataset_name}_val.xyz', index = ':', format = 'extxyz'),
                       desc = 'loading val  ', total = 1331):
            atoms_list_val.append(a)


        for a in tqdm(iread(f'{folder}/{dataset_name}_test.xyz', index = ':', format = 'extxyz'),
                       desc = 'loading test ', total = 1317):
            atoms_list_test.append(a)
        return atoms_list_train, atoms_list_val, atoms_list_test, index

    if dataset_name == 'nebBVSE122k':
        index = pd.read_csv(f'{folder}/{dataset_name}_index.csv')
        atoms_list_train, atoms_list_val, atoms_list_test = [], [], []
        for a in tqdm(iread(f'{folder}/{dataset_name}_train.xyz', index = ':', format = 'extxyz'),
                       desc = 'loading train', total = 96849):
            atoms_list_train.append(a)


        for a in tqdm(iread(f'{folder}/{dataset_name}_val.xyz', index = ':', format = 'extxyz'),
                       desc = 'loading val  ', total =  12405):
            atoms_list_val.append(a)


        for a in tqdm(iread(f'{folder}/{dataset_name}_test.xyz', index = ':', format = 'extxyz'),
                       desc = 'loading test ', total = 13167):
            atoms_list_test.append(a)
        return atoms_list_train, atoms_list_val, atoms_list_test, index



def load_index(dataset_name, folder):

    """
    Returns an index DataFrame

    Parameters
    ----------

    dataset_name: str
        dataset name, can be "BVEL13k", "nebBVSE122k", "nebDFT2k"

    folder: str
        path to the data folder
        

    Returns
    -------
    index: pandas DataFrame
        see examples

    Examples
    --------

    # nebDFT2k dataset
    >>> from litraj.data import load_index
    >>> df = load_index('nebDFT2k', 'pth/to/folder/data')

    """

    if dataset_name not in ['BVEL13k', 'nebBVSE122k', 'nebDFT2k']:
        raise ValueError(f'Wrong dataset name: {dataset_name}\nAllowed names: BVEL13k, nebBVSE122k, nebDFT2k')

    if dataset_name == 'BVEL13k':
        index = pd.read_csv(f'{folder}/{dataset_name}/{dataset_name}_index.csv')

    if dataset_name == 'nebBVSE122k':
        index = pd.read_csv(f'{folder}/{dataset_name}/{dataset_name}_index.csv')

    if dataset_name == 'nebDFT2k':
        index = pd.read_csv(f'{folder}/{dataset_name}/{dataset_name}_index.csv')

    return index



def _resource_path(relative_path):
    """Get absolute path to resource"""
    base_path = getattr(sys, '_MEIPASS', os.path.dirname(os.path.abspath(__file__)))
    path = os.path.join(base_path, relative_path)
    return path



def _download(url, filename):
    # credit: https://blasferna.com/articles/visualizing-download-progress-with-tqdm-in-python/
    try:
        response = requests.get(url, stream=True, allow_redirects=True)
        response.raise_for_status()  # Raises an exception for HTTP error status codes
        total_size = int(response.headers.get('content-length', 0))
        with open(filename, 'wb') as file, tqdm(
            desc=f'Downloading {filename}',
            total=total_size,
            unit='B',
            unit_scale=True,
            unit_divisor=1024,
        ) as bar:
            for data in response.iter_content(chunk_size=1024):
                size = file.write(data)
                bar.update(size)
    except RequestException as e:
        print(f"Error during download: {e}")
    except IOError as e:
        print(f"Error writing file: {e}")



def _unzip(pth, folder):
    with ZipFile(f'{pth}.zip', 'r') as zf:
        for member in tqdm(zf.infolist(), desc='Extracting'):
            try:
                zf.extract(member, folder)
            except Exception as err:
                raise err

