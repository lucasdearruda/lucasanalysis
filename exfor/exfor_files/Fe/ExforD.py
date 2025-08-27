import requests

def download_exfor(sub_id, internal_name):
    base_url = "https://www-nds.iaea.org/exfor/servlet/X4sGetSubent?reqx=98238&subID="
    url = base_url + sub_id.replace(".", "")  # Ensure URL formatting
    
    response = requests.get(url)
    
    if response.status_code == 200:
        filename = f"{internal_name}.txt"
        with open(filename, "w", encoding="utf-8") as file:
            file.write(response.text)
        print(f"File downloaded successfully as {filename}")
    else:
        print("Failed to download file. Check the ID and try again.")

if __name__ == "__main__":
    sub_id = input("Enter EXFOR Sub ID (e.g., 22718.002): ")
    internal_name = input("Enter internal name for the file: ")
    download_exfor(sub_id, internal_name)
