import os
from flask import Flask, flash, request, redirect, url_for, render_template, send_from_directory, send_file 
from werkzeug.utils import secure_filename
import config
from binding_site_search import search

ALLOWED_EXTENSIONS = {'xlsx', 'csv'}

path = os.path.dirname(os.path.abspath(__file__))
try:
    os.mkdir("uploads")
except FileExistsError:
    pass

upload_path = os.path.join(path, "uploads")
default_path = os.path.join(path, "defaults")
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

        uploaded_unbound = request.files['bound_file']
        uploaded_bound = request.files['unbound_file']
        uploaded_compound = request.files['compound_file']

        tolerance = request.form.get('tolerance')
        peak_height = request.form.get('peak_height')
        multi_protein = request.form.get('multiprotein')
        min_primaries = request.form.get('min_primaries')
        full_data = request.form.get('fulldata')

        default, data_dir = check_uploaded_files(uploaded_unbound.filename, uploaded_bound.filename, uploaded_compound.filename, analysis_complete)

        # file is submitted and file type is correct so now grab the file
        # to do - permission error when trying to save file to uploads folder
        if not default:
            download_uploaded_file(config.bound_filename, request, "bound_file")
            download_uploaded_file(config.unbound_filename, request, "unbound_file")
            download_uploaded_file(config.compounds_list_filename, request, "compound_file")

        # run external python module
        # create file paths to read the file
        bound_file_path = os.path.join(data_dir, config.bound_filename)
        unbound_file_path = os.path.join(data_dir, config.unbound_filename)
        compounds_file_path = os.path.join(data_dir, config.compounds_list_filename)

        binding_site_df = search(bound_file_path, unbound_file_path, compounds_file_path, \
            tolerance, peak_height, multi_protein, min_primaries, full_data)
        print(binding_site_df)

        # download
        outputs = os.path.join(path, download_path)
        output_csv = os.path.join(outputs, "BindingSites.csv")
        # output_html = os.path.join(outputs, "BindingSites.html")

        binding_site_df.to_csv(output_csv, index=False)
        # binding_site_df.to_html(output_html)
        analysis_complete = True

    return render_template("home.html", analysis_complete=analysis_complete)

def download_uploaded_file(filename, request, file_id):
    request.files[file_id].save(os.path.join(upload_path, filename))

@app.route('/download', methods=['GET', 'POST'])
def download(): 
    # Appending app path to upload folder path within app root folder
    outputs = os.path.join(path, download_path)
    print(outputs)
    # Returning file from appended path
    return send_file(os.path.join(outputs, 'BindingSites.csv'), as_attachment=True) ### NEED TO USE INCOGNITO ###

@app.after_request
def cache_control(response):
    response.headers['X-UA-Compatible'] = 'IE=Edge,chrome=1'
    response.headers['Cache-Control'] = 'public, max-age=0'
    return response

def check_uploaded_files(unbound, bound, compound, analysis_complete):
    '''
    Checks if files are uploaded correctly and returns True if no files uploaded
    ''' 
    print(unbound, bound, compound)
    default = False
    data_dir = upload_path

    if unbound == "" and bound == "" and compound == "":
        default = True
        data_dir = default_path
        flash("No files uploaded - default used")

    elif unbound == "" or bound == "" or compound == "":
        # nothing is uploaded
        flash("Please upload all required files")
        return render_template("home.html", analysis_complete=analysis_complete)

    else:
        # some file is entered
        # check it is the valid file type
        allowed_files_mask = allowed_file(unbound) and allowed_file(bound) and allowed_file(compound)
        if not allowed_files_mask:
            flash("Wrong File Type. Please only upload csv or xlsx files. ")
            return render_template("home.html", analysis_complete=analysis_complete)

    return default, data_dir

if __name__ == "__main__":
    app.run(debug=True)