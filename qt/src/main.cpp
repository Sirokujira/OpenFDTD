// main.cpp — OpenFDTD-X application entry
#include <QApplication>
#include <QCommandLineParser>
#include <QFile>
#include <QSettings>
#include <QStyleFactory>
#include <QTimer>

#include "MainWindow.h"
#include "I18n.h"

int main(int argc, char *argv[])
{
    QApplication app(argc, argv);
    QApplication::setOrganizationName("OpenFDTD");
    QApplication::setApplicationName("OpenFDTD-X");
    QApplication::setApplicationVersion("1.0.0");

    // Force Fusion style on all platforms so the look matches the mock
    // (Windows Vista / macOS native styles diverge too much).
    QApplication::setStyle(QStyleFactory::create("Fusion"));

    QCommandLineParser cli;
    cli.setApplicationDescription("OpenFDTD-X — Multi-Domain FDTD GUI");
    cli.addHelpOption();
    cli.addVersionOption();
    cli.addPositionalArgument("file", "Open a .ofd project file", "[file]");
    QCommandLineOption langOpt({ "l", "lang" }, "UI language (ja|en|both)", "lang");
    cli.addOption(langOpt);
    QCommandLineOption domainOpt({ "d", "domain" },
        "Start in domain (em|optical|acoustic|underwater)", "domain");
    cli.addOption(domainOpt);
    QCommandLineOption shotOpt("screenshot",
        "Save a window screenshot to <path> and exit (for CI)", "path");
    cli.addOption(shotOpt);
    cli.process(app);

    // i18n: CLI option > saved setting > ja
    const QString lang = cli.isSet(langOpt)
        ? cli.value(langOpt)
        : QSettings().value("ui/language", "ja").toString();
    ofd::I18n::instance().setLanguage(lang);

    // stylesheet matching the HTML mock's Qt Fusion look
    QFile qss(":/styles/openfdtd.qss");
    if (qss.open(QIODevice::ReadOnly | QIODevice::Text))
        app.setStyleSheet(QString::fromUtf8(qss.readAll()));

    ofd::MainWindow w;
    w.show();

    const QStringList args = cli.positionalArguments();
    if (!args.isEmpty())
        w.openProject(args.first());

    if (cli.isSet(domainOpt))
        w.setDomain(ofd::domainFromKey(cli.value(domainOpt)));

    if (cli.isSet(shotOpt)) {
        const QString path = cli.value(shotOpt);
        QTimer::singleShot(800, &w, [&w, path] {
            w.grab().save(path);
            QApplication::quit();
        });
    }

    return app.exec();
}
