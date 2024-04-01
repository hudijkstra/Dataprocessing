def check_directory(directory, extension):
    for filename in os.listdir(directory):
        print(filename)
        if filename.endswith(extension):
            return os.path.join(directory, filename)
    return None
