import os
from flask import Flask 
from flask import render_template,request,redirect,url_for,flash
from werkzeug.utils import secure_filename
def create_app(test_config=None):
    # create and config the app
    app = Flask(__name__, instance_relative_config=True)
    app.config.from_mapping(
        SECRET_KEY='dev',
        DATABASE=os.path.join(app.instance_path, 'flaskr.sqlite'),
        UPLOAD_FOLDER = os.getcwd()+'/MVP/upload_files'
    )
    #os.makedirs(os.path.join(app.instance_path,'upload_files'), exist_ok=True)

    if test_config is None:
        # load the instance config, if it exists, when not testing
        app.config.from_pyfile('config.py', silent=True)
    else:
        # load the test config if passed in
        app.config.from_mapping(test_config)

    # ensure the instance folder exists
    try:
        os.makedirs(app.instance_path)
    except OSError:
        pass

    # a simple page that says hello
    @app.route('/hello')
    def hello():
        return 'Hello, World!'
    @app.route('/',methods=['GET','POST'])
    def upload_file():
        if request.method == 'POST':
            # check if the post request has the file part
            if 'file' not in request.files:
                flash('No file part')
                return redirect(request.url)
            file = request.files['file']
            # if user does not select file, browser also
            # submit an empty part without filename
            if file.filename == '':
                flash('No selected file')
                return redirect(request.url)
            if file :
                filename = secure_filename(file.filename)
                file.save(os.path.join(app.config['UPLOAD_FOLDER'],filename))
                #return app.config['UPLOAD_FOLDER']+" file uploaded "+filename
                #return redirect(url_for('uploaded_file',
                                    #filename=filename))
        return render_template('base.html')
        #'''<!doctype html>
        #<title>Upload new File</title>
        
        #</form>
        #'''
        



    from . import graph
    app.static_file='static'
    app.register_blueprint(graph.bp)
    return app
