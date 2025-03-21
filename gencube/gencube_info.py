from pathlib import Path      # For handling filesystem paths in a platform-independent way

# Custom functions
from .utils import (
    get_entrez_info,
    delete_file,
)

## gencube info ------------------------------
def info (
    info_save
):

    if not info_save:
        home_dir = Path.home()
        config_file = home_dir / '.gencube_entrez_info'

        print(
            'The submitted email and API key are as follows.\n'
            '\n'
            f'Path: {config_file}'
            )

        # Show the previously submitted information
        with open(config_file, 'r') as in_f:
            for line in in_f:
                if '#' not in line:
                    print(line.replace('\n', ''))
            print('')

        while True:
            question = input("Would you like to modify them? (type 'yes or 'no'): ").strip()

            if question.lower() in ['yes', 'y']:
                delete_file(config_file)
                get_entrez_info()
                break

            elif question.lower() in ['no', 'n']:
                break




