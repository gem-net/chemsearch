import os
import logging

basedir = os.path.abspath(os.path.dirname(__file__))
_logger = logging.getLogger(__name__)


class Config:
    SECRET_KEY = os.environ.get('SECRET_KEY') or 'hard to guess string'

    USE_AUTH = os.environ.get('USE_AUTH', 'off').lower() not in {'off', 'false', '0'}
    OAUTH_CREDENTIALS = {
        'google': {
            'id': os.environ.get('GOOGLE_CLIENT_ID'),
            'secret': os.environ.get('GOOGLE_SECRET')
        }
        }

    MOLECULES_PER_PAGE = os.environ.get('MOLECULES_PER_PAGE', 15)

    SQLALCHEMY_DATABASE_URI = \
        'sqlite:///' + os.path.join(basedir, 'db.sqlite')
    SQLALCHEMY_TRACK_MODIFICATIONS = False
    # SESSION_TYPE = 'filesystem'
    # SESSION_FILE_DIR = os.environ.get('SESSION_FILE_DIR', os.getcwd())
    # _logger.info(f"{SESSION_FILE_DIR=}")
    MAX_CONTENT_LENGTH = 1024 * 1024  # 1MB request size limit, for uploads
    MAIL_SERVER = os.environ.get('MAIL_SERVER', 'smtp.googlemail.com')
    MAIL_PORT = int(os.environ.get('MAIL_PORT', '587'))
    MAIL_USE_TLS = os.environ.get('MAIL_USE_TLS', 'true').lower() in \
        ['true', 'on', '1']
    MAIL_USERNAME = os.environ.get('MAIL_USERNAME')
    MAIL_PASSWORD = os.environ.get('MAIL_PASSWORD')
    MAIL_SUBJECT_PREFIX = '[Molecules]'
    MAIL_SENDER = os.environ.get('MAIL_SENDER')
    MAIL_ADMIN = os.environ.get('MAIL_ADMIN')

    SHARED_DRIVE_ID = os.environ.get('SHARED_DRIVE_ID')

    LOCAL_DB_PATH = os.environ.get('LOCAL_DB_PATH')

    @staticmethod
    def init_app(app):
        pass


class DevelopmentConfig(Config):
    DEBUG = True


class TestingConfig(Config):
    TESTING = True


class ProductionConfig(Config):
    @classmethod
    def init_app(cls, app):
        print("Running the production config init code.")
        Config.init_app(app)

        # email errors to the administrators
        import logging
        from logging.handlers import SMTPHandler
        credentials = None
        secure = None
        if getattr(cls, 'MAIL_USERNAME', None) is not None:
            credentials = (cls.MAIL_USERNAME, cls.MAIL_PASSWORD)
            if getattr(cls, 'MAIL_USE_TLS', None):
                secure = ()
        # print(f"{cls.MAIL_SERVER=}, {cls.MAIL_ADMIN=}, {cls.MAIL_SERVER}")
        mail_handler = SMTPHandler(
            mailhost=(cls.MAIL_SERVER, cls.MAIL_PORT),
            fromaddr=cls.MAIL_SENDER,
            toaddrs=[cls.MAIL_ADMIN],
            subject=cls.MAIL_SUBJECT_PREFIX + ' Application Error',
            credentials=credentials,
            secure=secure)
        mail_handler.setLevel(logging.ERROR)
        app.logger.addHandler(mail_handler)


config = {
    'development': DevelopmentConfig,
    'testing': TestingConfig,
    'production': ProductionConfig,

    'default': DevelopmentConfig
}
