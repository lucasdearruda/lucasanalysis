import requests

def get_repos_size():
    token = input("Enter your GitHub token: ").strip()
    headers = {'Authorization': f'token {token}'}

    total_size_kb = 0
    page = 1
    while True:
        url = f'https://api.github.com/user/repos?per_page=100&page={page}'
        r = requests.get(url, headers=headers)

        if r.status_code != 200:
            print(f"Error: {r.status_code} - {r.text}")
            break

        repos = r.json()
        if not repos:
            break

        for repo in repos:
            size_kb = repo.get('size', 0)
            total_size_kb += size_kb
            print(f"{repo['name']}: {size_kb} KB")

        page += 1

    print(f"Total size: {total_size_kb / 1024:.2f} MB")

if __name__ == "__main__":
    get_repos_size()
