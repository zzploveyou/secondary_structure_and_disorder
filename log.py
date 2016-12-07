# coding:utf-8
'''
Author       :  ZZP
Mail         :  zhangzhaopeng@mail.nankai.edu.cn
Created Time :  2016-11-16 19:12:09

Description  :  log sth with tags and colors using ASSIC.

'''

def gettime(format="%m-%d %H:%M:%S"):
    import datetime
    return datetime.datetime.strftime(datetime.datetime.now(), format)

class Terminal_log:
    """
    log into screen with tags and colors.

    Members
    -------
    debug(self, msg)
        debug message.
    info(self, msg)
        info message.
    warn(self, msg)
        warn message.
    error(self, msg)
        error message.
    fatal(self, msg)
        fatal error message.
    done(self, msg)
        done message.
    """

    def __init__(self, logfilename=None, tag=True, color=True, time=True, timeformat="%H:%M:%S", brief=False):
        """
        constructor of Terminal_log.

        Parameters
        ----------
        logfilename: str
            default: None
            log into file at the same time of log into screen if it's true.
        tag: bool
            default: True
            log without tags if it's false.
        color: bool
            default: True
            log without colors if it's false.
        time: bool
            default: True
            log with time if it's true.
        timeformat: str
            default: "%m-%d %H:%M:%S"
            time format for log time.
        brief: bool
            default: False
            don't print debug log if it's true.
        """
        if color:
            self.__ENDC__ = "\033[0m"
            self.__INFO__ = "\033[32m"
            self.__WARN__ = "\033[33m"
            self.__ERROR__ = "\033[35m"
            self.__FATAL__ = "\033[31m"
            self.__DONE__ = "\033[36m"
        else:
            self.__ENDC__ = ""
            self.__INFO__ = ""
            self.__WARN__ = ""
            self.__ERROR__ = ""
            self.__FATAL__ = ""
            self.__DONE__ = ""
        
        self._logfilename = logfilename
        self._tag = tag
        self.time = time
        self.timeformat = timeformat
        self._brief = brief
        
        if self._logfilename:
            fw = open(self._logfilename, 'a')
            fw.write('############### %s ###############\n' %gettime(format=self.timeformat))
            fw.close()

    def __logfw__(self, msg):
        """
        write to self._logfilename

        Parameters
        ----------
        msg: str
            message that would be written to self._logfilename
        """
        fw = open(self._logfilename, 'a')
        fw.write(msg + "\n")
        fw.close()

    def debug(self, msg):
        """
        log debug message.

        Parameters
        ----------
        msg: str
            print debug message.
            if self._brief == True, don't print.

        Examples
        --------
        >>> mylog = Terminal_log()
        >>> mylog.debug("this is debug")
        [d] this is debug
        """
        if not self._brief:
            if self.time:
                msg = "%s %s" %(gettime(format=self.timeformat), msg)
            msg = "[d] " + msg
            print msg
            if self._logfilename:
                self.__logfw__(msg)

    def info(self, msg):
        """
        log info message.

        Parameters
        ----------
        msg: str
            print info message.

        Examples
        --------
        >>> mylog = Terminal_log()
        >>> mylog.info("this is info")
        \033[32m[+] this is info\033[0m
        """
        if self.time:
            msg = "%s %s" %(gettime(format=self.timeformat), msg)
        if self._tag:
            msg = "[+] " + msg
        result = self.__INFO__ + msg + self.__ENDC__
        print(result)
        if self._logfilename:
            self.__logfw__(msg)

    def warn(self, msg):
        """
        log warn message.

        Parameters
        ----------
        msg: str
            print warn message.

        Examples
        --------
        >>> mylog = Terminal_log()
        >>> mylog.warn("this is warn")
        \033[33m[!] this is warn\033[0m
        """
        if self.time:
            msg = "%s %s" %(gettime(format=self.timeformat), msg)
        if self._tag:
            msg = "[!] " + msg
        result = self.__WARN__ + msg + self.__ENDC__
        print(result)
        if self._logfilename:
            self.__logfw__(msg)

    def error(self, msg):
        """
        log error message.

        Parameters
        ----------
        msg: str
            print error message.

        Examples
        --------
        >>> mylog = Terminal_log()
        >>> mylog.error("this is error")
        \033[31m[-] this is error\033[0m
        """
        if self.time:
            msg = "%s %s" %(gettime(format=self.timeformat), msg)
        if self._tag:
            msg = "[-] " + msg
        result = self.__ERROR__ + msg + self.__ENDC__
        print(result)
        if self._logfilename:
            self.__logfw__(msg)

    def fatal(self, msg):
        """
        log fatal error message.

        Parameters
        ----------
        msg: str
            print fatal error message.

        Examples
        --------
        >>> mylog = Terminal_log()
        >>> mylog.fatal("this is fatal error")
        \033[35m[×] this is fatal error\033[0m
        """
        if self.time:
            msg = "%s %s" %(gettime(format=self.timeformat), msg)
        if self._tag:
            msg = "[×] " + msg
        result = self.__FATAL__ + msg + self.__ENDC__
        print(result)
        if self._logfilename:
            self.__logfw__(msg)

    def done(self, msg):
        """
        log done message.

        Parameters
        ----------
        msg: str
            print done message.

        Examples
        --------
        >>> mylog = Terminal_log()
        >>> mylog.done("done.")
        \033[35m[✓] done.\033[0m
        """
        if self.time:
            msg = "%s %s" %(gettime(format=self.timeformat), msg)
        if self._tag:
            msg = "[✓] " + msg
        result = self.__DONE__ + msg + self.__ENDC__
        print(result)
        if self._logfilename:
            self.__logfw__(msg)

if __name__ == '__main__':
    """
    test Terminal_log
    """
    m = Terminal_log(color=True, timeformat="%H:%M:%S")
    m.debug("this is debug")
    m.info("this is info")
    m.warn("this is warn")
    m.error("this is error")
    m.fatal("this is fatal")
    m.done("sth is done.")
