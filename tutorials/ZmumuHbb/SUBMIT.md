# SubMIT
MIT provides computing resources specifically tailored for physics research, including High Energy Physics (HEP). The main computing cluster is called SubMIT:  [https://submit.mit.edu/](https://submit.mit.edu/)

SubMIT supports a variety of projects, including studies for the Future Circular Collider (FCC).

In this section, we’ll walk through the basics of: creating an account, logging into SubMIT and using its computing resources for analysis.



## MIT account
The first step is to obtain an MIT account to access SubMIT and related services. For non-MIT or external users, the following information is required:

- First and last name  
- Date of birth  
- Institution  
- Email address  
- Phone number  

This information should be sent to the designated MIT contact person. Once submitted, your account will be created and access credentials will be provided.


## Access to the SubMIT cluster
To access the cluster, follow the steps outlined in the  ["Getting Started" section](https://submit.mit.edu/submit-users-guide/starting.html) of the  [SubMIT User's Guide](https://submit.mit.edu/submit-users-guide/index.html). 
Authentication requires the use of **SSH keys**, which are explained in the guide.

Once your keys are set up, you can log in to the SubMIT cluster using the following command in your terminal:

```sh
ssh <your_username>@submit.mit.edu
```


## JupyterHub
You can also access subMIT via JupyterHub, which provides a web-based interface to the cluster. You can access JupyterHub [https://submit.mit.edu/jupyter](https://submit.mit.edu/jupyter). Documentation for this can be found [in the subMIT User's Guide](https://submit.mit.edu/submit-users-guide/access.html#jupyterhub).


## VS Code
We suggest to use Visual Studio Code (VS Code) as a text editor. It is a powerful and user-friendly editor that is widely used in the scientific community. You can download it [here](https://code.visualstudio.com/). It also has a built in command line, so that you can easily use your code all in one application. You may also do all the following steps there if you wish.

Instructions to set up VS Code on the subMIT cluster can be found on the [subMIT User's Guide](https://submit.mit.edu/submit-users-guide/access.html#vscode).


## Personal Website on subMIT
As documented in the [subMIT User's Guide](https://submit.mit.edu/submit-users-guide/starting.html#creating-a-personal-webpage), you can create a directory called `public_html` in your home directory on the subMIT cluster. Any files you put in this directory will be accessible on the web at `http://submit08.mit.edu/~<username>/`. This is a great way to share your scripts, plots, etc. with others.

You can add your own .php files to your `public_html` directory to edit the style of your webpage.
