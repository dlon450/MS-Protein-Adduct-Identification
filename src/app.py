import os
from flask import Flask, flash, request, redirect, url_for, render_template, send_from_directory, send_file 
from werkzeug.utils import secure_filename
import config
from Binding_Site_Search import search

ALLOWED_EXTENSIONS = {'xlsx', 'csv'}

path = os.path.dirname(os.path.abspath(__file__))
try:
    os.mkdir("uploads")
except FileExistsError:
    pass

upload_path = os.path.join(path, "uploads")
download_path = os.path.join(path, "outputs")

def allowed_file(filename):
    '''
    Check file is of valid type (.xlsx or .csv)
    '''
    return "." in filename and filename.rsplit(".", 1)[1].lower() in ALLOWED_EXTENSIONS

app = Flask(__name__)
app.secret_key = b'_5#y2L"F4Q8z\n\xec]/'

# home page
@app.route("/")
def home():
    analysis_complete = False
    return render_template("home.html",analysis_complete=analysis_complete)
    
@app.route("/upload", methods=["GET", "POST"])
def upload():
    analysis_complete = False

    if request.method == "POST":
        # grab filename
        unbound_filename = config.unbound_filename
        bound_filename = config.bound_filename
        compound_file = config.compounds_list_filename

        if unbound_filename == "" or bound_filename == "" or compound_file == "":
            # nothing is uploaded
            flash("Please upload all required files")
            return render_template("home.html", analysis_complete=analysis_complete)

        # some file is entered
        # check it is the valid file type
        allowed_files_mask = allowed_file(unbound_filename) and allowed_file(bound_filename) and allowed_file(compound_file)
        if not allowed_files_mask:
            flash("Wrong File Type. Please only upload xlsx files. ")
            return render_template("home.html", analysis_complete=analysis_complete)

        # file is submitted and file type is correct so now grab the file
        # to do - permission error when trying to save file to uploads folder
        download_uploaded_file(bound_filename, request, "bound_file")
        download_uploaded_file(unbound_filename, request, "unbound_file")
        download_uploaded_file(compound_file, request, "compound_file")

        print('Downloaded!')

        # run external python module
        # create file paths to read the file
        bound_file_path = os.path.join(upload_path, bound_filename)
        unbound_file_path = os.path.join(upload_path, unbound_filename)
        compounds_file_path = os.path.join(upload_path, compound_file)

        binding_site_df = search(bound_file_path, unbound_file_path, compounds_file_path)

        print('Searching complete!')

        # download
        outputs = os.path.join(path, download_path)
        outputs = os.path.join(outputs, "BindingSites.xlsx")
        binding_site_df.to_excel(outputs, index=False)
        analysis_complete = True

    return render_template("home.html", analysis_complete=analysis_complete)

def download_uploaded_file(filename, request, file_id):
    request.files[file_id].save(os.path.join(upload_path, filename))

@app.route('/download/<filename>', methods=['GET', 'POST'])
def download(filename): 
    # Appending app path to upload folder path within app root folder
    outputs = os.path.join(path, download_path)
    # Returning file from appended path
    return send_from_directory(directory=outputs, filename=filename,as_attachment=True)
    
    
if __name__ == "__main__":
    app.run(debug=True)