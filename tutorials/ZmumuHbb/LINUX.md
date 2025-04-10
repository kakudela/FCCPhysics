
# The command line

Login into the subMIT cluster  using this command in the terminal:
```sh
ssh <your_username>@submit.mit.edu
```
If you see this at the beginning of your line
```sh
[<your_username>@submit{a number} ~]$
```
you're good to go and successfully logged in!

You will do your work here from now on, rather than on your laptop.

We are ready to use the terminal. Finally, you can be just like Mr. Robot, and impress all your family and friends. Covering the full of use of the command line is far outside the scope of this tutorial, but we will cover the basics. A more comprehensive guide can be found [here](https://ubuntu.com/tutorials/command-line-for-beginners#3-opening-a-terminal).

Basic commands:

- `ls`: List the contents of the current directory.
- `cd`: Change the current directory.
- `pwd`: Print the current directory.
- `mkdir`: Create a new directory.
- `rm`: Remove a file (warning: there is no Trash, it will be forever deleted).
- `cp`: Copy a file.
- `mv`: Move a file.
- `cat`: Print the contents of a file.
- `vim`: Open a file in the vim text editor.
- `history`: Show a list of previous commands.

Basic shortcuts:
- `tab`: Autocomplete a command or file name.
- `up arrow`: Scroll through previous commands.
- `ctrl + c`: Stop a running command.
- `ctrl + d`: Exit the terminal.
- `ctrl + l`: Clear the terminal.
- `ctrl + r`: Search through previous commands.

> *Exercise*: In your home space, create a new directory called `fcc-ee`, and navigate to it. Print the full path of your current directory.


# Git

Git is a version control system that allows you to track changes in your code as well. It is widely used in the scientific community and is essential for collaborative work. We will use Git to manage our code. In particular, we will use GitHub, a platform that hosts Git repositories.

> *Exercise*: Create a GitHub account. Setup your GitHub keys on subMIT. Fork this repository (https://github.com/mit-fcc/tutorials). Navigate to the directory you created earlier, `fcc-ee`, and clone the repository. Edit the README file by adding your name and the project you are working on, and push the changes to your fork. Navigate to the repository on GitHub's website and check that the changes are there. You will need the following commands for this: `git clone`, `git add`, `git commit`, and `git push` and you can have detailed introduction on the HSF websites [here](https://hsf-training.github.io/analysis-essentials/git/README.html).
