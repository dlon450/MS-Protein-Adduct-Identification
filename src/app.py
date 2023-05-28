import os
from flask import Flask, flash, request, render_template, send_file 
# from werkzeug.utils import secure_filename
import config
from binding_site_search import search

ALLOWED_EXTENSIONS = {'xlsx', 'csv'}

path = os.path.dirname(os.path.abspath(__file__))

upload_path = os.path.join(os.path.dirname(path), "uploads")
default_path = os.path.join(os.path.dirname(path), "example")
download_path = os.path.join(os.path.dirname(path), "outputs")

try:
    os.mkdir(upload_path)
    os.mkdir(download_path)
except FileExistsError:
    pass

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
        uploaded_compound = request.files['compound_file']
        uploaded_adducts = request.files['adducts_file']
        if uploaded_adducts.filename == "":
            adducts_file_path = os.path.join(default_path, config.adducts_filename)
        else:
            download_uploaded_file(config.adducts_filename, request, "adducts_file")
            adducts_file_path = os.path.join(upload_path, config.adducts_filename)

        tolerance = request.form.get('tolerance')
        peak_height = request.form.get('peak_height')
        multi_protein = 'on'
        min_primaries = request.form.get('min_primaries')
        max_primaries = request.form.get('max_primaries')
        max_adducts = request.form.get('max_adducts')
        valence = request.form.get('valence')
        min_dist_between_peaks = request.form.get('min_dist_between_peaks')
        calibrate = request.form.get('calibrate')
        manual_calibration = request.form.get('manual_calibration_amount')
        only_best = 'off'
        return_all_peaks = request.form.get('return_all_peaks')
        isotope_pattern_method = request.form.get('isotope_pattern_method')

        default, data_dir = check_uploaded_files(uploaded_unbound.filename, uploaded_compound.filename, analysis_complete)

        # file is submitted and file type is correct so now grab the file
        # to do - permission error when trying to save file to uploads folder
        if not default:
            download_uploaded_file(config.bound_filename, request, "bound_file")
            download_uploaded_file(config.compounds_list_filename, request, "compound_file")

        # run external python module
        # create file paths to read the file
        bound_file_path = os.path.join(data_dir, config.bound_filename)
        compounds_file_path = os.path.join(data_dir, config.compounds_list_filename)

        binding_site_df = search(bound_file_path, compounds_file_path, adducts_file_path, tolerance, \
            peak_height, multi_protein, min_primaries, max_primaries, max_adducts, valence, only_best, \
                min_dist_between_peaks, calibrate, manual_calibration)
        print(binding_site_df)

        # download
        outputs = os.path.join(path, download_path)
        output_csv = os.path.join(outputs, "BindingSites.csv")

        binding_site_df.to_csv(output_csv, index=False)
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

def check_uploaded_files(bound, compound, analysis_complete):
    '''
    Checks if files are uploaded correctly and returns True if no files uploaded
    ''' 
    print(bound, compound)
    default = False
    data_dir = upload_path

    if bound == "" and compound == "":
        default = True
        data_dir = default_path
        flash("No files uploaded - default used")

    elif bound == "" or compound == "":
        # nothing is uploaded
        flash("Please upload all required files")
        return render_template("home.html", analysis_complete=analysis_complete)

    else:
        # some file is entered
        # check it is the valid file type
        allowed_files_mask = allowed_file(bound) and allowed_file(compound)
        if not allowed_files_mask:
            flash("Wrong File Type. Please only upload csv or xlsx files. ")
            return render_template("home.html", analysis_complete=analysis_complete)

    return default, data_dir

if __name__ == "__main__":
    # app.run(debug=True)
    # app.run()
    from waitress import serve
    serve(app, port=8080)