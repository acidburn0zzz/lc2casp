#include "translator.hh"
#include "printer.hh"
#include <potassco/convert.h>
#include <potassco/aspif.h>
#include <program_opts/application.h>
#include <program_opts/typed_value.h>
#include <fstream>
#include <iostream>
#include <cctype>

using namespace ProgramOptions;

class LpConvert : public ProgramOptions::Application {
public:
    virtual const char* getName()       const { return "lc2casp"; }
    virtual const char* getVersion()    const { return "1.0.0"; }
    virtual PosOption   getPositional() const { return &positional; }
    virtual const char* getUsage()      const {
        return
            "[options] [<file>]\n"
            "Translate program in <file> or standard input";
    }
    virtual void initOptions(OptionContext& root);
    virtual void validateOptions(const OptionContext&, const ParsedOptions&, const ParsedValues&) {}
    virtual void setup() {}
    virtual void run();
private:
    static bool positional(const std::string&, std::string& optOut) {
        optOut = "input";
        return true;
    }
    static int error(int line, const char* what) {
        fprintf(stderr, "*** ERROR: In line %d: %s\n", line, what);
        static_cast<LpConvert*>(Application::getInstance())->exit(EXIT_FAILURE);
        return EXIT_FAILURE;
    }
    std::string input_;
    std::string output_;
    std::pair<int, int> bound_ = {std::numeric_limits<int>::min(), std::numeric_limits<int>::max()};
    bool text_ = false;
};

void LpConvert::initOptions(OptionContext& root) {
    OptionGroup convert("Options");
    convert.addOptions()
        ("input,i,@2", storeTo(input_), "Input file")
        ("text,t", storeTo(text_)->flag(), "Do not translate but print the input in something more readable")
        ("bounds,b", storeTo(bound_), "Pair of values limiting the minimum and maximum value for integer variables")
        ("output,o", storeTo(output_)->arg("<file>"), "Write output to <file> (default: stdout)")
    ;
    root.add(convert);
}
void LpConvert::run() {
    std::ifstream iFile;
    std::ofstream oFile;
    if (!input_.empty() && input_ != "-") {
        iFile.open(input_.c_str());
        if (!iFile.is_open()) { throw std::runtime_error("Could not open input file!"); }
    }
    if (!output_.empty() && output_ != "-") {
        if (input_ == output_) { throw std::runtime_error("Input and output must be different!"); }
        oFile.open(output_.c_str());
        if (!oFile.is_open()) { throw std::runtime_error("Could not open output file!"); }
    }
    std::istream& in = iFile.is_open() ? iFile : std::cin;
    std::ostream& os = oFile.is_open() ? oFile : std::cout;
    if (in.peek() == 'a') {
        if (text_) {
            Printer writer(os);
            Potassco::AspifInput reader(writer);
            readProgram(in, reader, error);
        }
        else {
            FoundedOutput writer(os, bound_.first, bound_.second);
            Potassco::AspifInput reader(writer);
            readProgram(in, reader, error);
        }
    }
    else {
        throw std::runtime_error("Unrecognized input format!");
    }
    iFile.close();
    oFile.close();
}

int main(int argc, char** argv) {
    LpConvert app;
    return app.main(argc, argv);
}
