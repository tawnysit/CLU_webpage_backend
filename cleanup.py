import os
from datetime import date, datetime, timedelta

def del_old_files(dir_path, save_days=14):
    """
    Function to delete files older than save_days in dir_path.
    """
    date_to_delete = date.today() - timedelta(days=save_days)
    
    for file in os.listdir(dir_path):
        full_path = os.path.join(dir_path, file)
        modtime = date.fromtimestamp(os.stat(full_path).st_mtime)
        
        if modtime < date_to_delete:
            if os.path.isdir(full_path): # check if directory/folder and use rmdir if yes
                try:
                    os.rmdir(full_path)
                except OSError:
                    print(f"Unable to remove dir: {full_path}")
                    
            if os.path.exists(full_path):
                try:
                    os.remove(full_path)
                except OSError:
                    print(f"Unable to remove file: {full_path}")
                    
def del_old_files_by_date(dir_path, save_date):
    """
    Function to delete files last modified before save_date in dir_path.
    save_date is a string in ISO format (YYYY-MM-DD)
    """
    date_to_delete = date.fromisoformat('2021-06-30')
    
    for file in os.listdir(dir_path):
        full_path = os.path.join(dir_path, file)
        modtime = date.fromtimestamp(os.stat(full_path).st_mtime)
        
        if modtime < date_to_delete:
            if os.path.isdir(full_path): # check if directory/folder and use rmdir if yes
                try:
                    os.rmdir(full_path)
                except OSError:
                    print(f"Unable to remove dir: {full_path}")
                    
            if os.path.exists(full_path):
                try:
                    os.remove(full_path)
                except OSError:
                    print(f"Unable to remove file: {full_path}")